(ns genome.annotate
  (require [clojure.java.io   :as io ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]
           [incanter.core     :as i  ]
           [genome.annotate   :as ga ]
           [incanter.io       :as ii ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRANGING A GFF3 FILE INTO A DATASET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn pairs [file]
  "parses through the ID line of gff3 file and creates a map "
  (->> file
         (re-seq  #"[^;]*")
         (map #(re-seq #"[^=]*" %))
         (map (fn [x] (filter #(not= "" %) x)))
         (filter #(= 2 (count %)))
         flatten
         (apply hash-map)))

(defn only-existing [id file]
  "puts a value if there is one if not a -"
  (let [value ((pairs file) id)]
    (if (= nil value)
      "-"
      value)))

(defn  add-id [id file]
  "here you choose what attribute you would like to parse"
  (->> file
       (i/add-derived-column
        (keyword id)
        [:attributes]
        #(only-existing id %))))

(defn get-list [file attribution]
  "create a file of the relevant data attribution :type/:gene"
  (->> (i/sel file :cols [attribution :pt1 :pt2])
       i/to-vect))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATING AN EMPTY REFFERENCE DATASET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-col-name [file attribute]
  "gets a vector of distinct col-names"
  (let [pos (->> file
                 (i/$ attribute)
                 distinct
                 (map #(str % "+"))
                 (map keyword)
                 vec)
        neg (->> file
                 (i/$ attribute)
                 distinct
                 (map #(str % "-"))
                 (map keyword)
                 vec)]
    (concat pos neg [:loc :gene+ :gene-])))

(defn empty-ref-db [file attribute size]
  "creates an empty ref dataset with col-names"
  (let [header (get-col-name file attribute)]
    (->> (i/dataset header (->> (repeat "-")
                                (take (count header))
                                vec
                                repeat
                                (take size)
                                vec))
         (i/add-column
          :loc
          (range 1 (inc size))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FILLING THE REFFERENCE DATASET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn add-data [new-type old-type loc pt1 pt2]
  "returns colname if in range"
  (if (and (>= loc pt1) (<= loc pt2))
    (name new-type)
    old-type))

(defn add-annotations [list posneg file]
  "adds annotations to empty ref dataset"
  (let [type (keyword (str (first list) posneg))
        pt1 (second list)
        pt2 (last list)]
    (->> file
         (i/add-derived-column
          type
          [type :loc]
          #(add-data type %1 %2 pt1 pt2)))))

(defn add-annotations-gene [list nontype file]
  "adds annotations to empty ref dataset"
  (let [type (keyword (first list))
        pt1 (second list)
        pt2 (last list)]
    (->> file
         (i/add-derived-column
          nontype
          [nontype :loc]
          #(add-data type %1 %2 pt1 pt2)))))

(defn creat-one [list posneg file]
  (let [new-file (add-annotations (first list) posneg file)
        new-list (rest list)]
    (if (< 1 (count list))
      (recur new-list posneg new-file)
      file)))

(defn creat-one-gene [list nontype file]
  (let [new-file (add-annotations-gene (first list) nontype file)
        new-list (rest list)]
    (if (< 1 (count list))
      (recur new-list nontype new-file)
      file)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;OPPORATION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn gff3>dataset [file]
  (let [reference  (ii/read-dataset file :delim \tab)
        renamed    (i/rename-cols {:col0 :seqid    :col1 :source
                                   :col2 :type     :col3 :pt1
                                   :col4 :pt2      :col5 :score
                                   :col6 :strand   :col7 :phase
                                   :col8 :attributes} reference)
        attributed (add-id "gene" renamed)
        scrubed    (->> attributed (i/$ [:seqid :gene :type :strand :pt1 :pt2 ]))
        scrubed+   (i/$where {:strand {:$eq "+"}} scrubed)
        scrubed-   (i/$where {:strand {:$eq "-"}} scrubed)
        rescrubed  (i/$where {:type {:$eq "gene"}} scrubed)
        rescrubed+ (i/$where {:strand {:$eq "+"}} rescrubed)
        rescrubed- (i/$where {:strand {:$eq "-"}} rescrubed)]
    
    (->> (empty-ref-db scrubed :type  (apply max (i/$ :pt2 scrubed)))
         (creat-one      (get-list scrubed+   :type)  "+")
         (creat-one      (get-list scrubed-   :type)  "-")
         (creat-one-gene (get-list rescrubed+ :gene) :gene+)
         (creat-one-gene (get-list rescrubed- :gene) :gene-))))

(defn gff3>csv [file_in file_out]
  "save annotation as csv file"
  (let [incfile (gff3>dataset file_in)]
    (with-open [f-out (io/writer file_out)]
      (csv/write-csv f-out [(map name (i/col-names incfile))])
      (csv/write-csv f-out (i/to-list incfile)))))

;("/home/yosh/datafiles/genes/merlin.gff3" "/home/yosh/datafiles/genes/merlin.inc")

(defn get-annotation [input_file]
  "creates an annotation"
  (let [ann_cols (ga/gff3>dataset input_file)]
    (->> (i/$ (vec (sort (distinct  (i/col-names ann_cols)))) ann_cols)
         (i/$ [:loc :gene+ :gene- :CDS+ :CDS- :exon- :exon+]))))
  (def m-get-annotation (memoize get-annotation))

(defn annotate [gff3 incfile]
  "incorporates annotation into the inc file"
  (let [annotation (m-get-annotation gff3)]
    (i/$join [:loc :ref-loc] annotation incfile)))
