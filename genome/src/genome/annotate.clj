(ns genome.annotate
  (require [clojure.java.io :as io]
           [clojure.string :as s]
           [incanter.core :as i]
           [incanter.stats :as st]
           [incanter.io :as ii ]
           [clojure.xml :as xml]))

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
    (println pos)
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

(defn add-annotations [file list posneg]
  "adds annotations to empty ref dataset"
  (let [type (keyword (str (first list) posneg))
        pt1 (second list)
        pt2 (last list)]
    (->> file
         (i/add-derived-column
          type
          [type :loc]
          #(add-data type %1 %2 pt1 pt2)))))

(defn add-annotations-gene [file list nontype]
  "adds annotations to empty ref dataset"
  (let [type (keyword (first list))
        pt1 (second list)
        pt2 (last list)]
    (->> file
         (i/add-derived-column
          nontype
          [nontype :loc]
          #(add-data type %1 %2 pt1 pt2)))))


(defn creat-one [file list posneg]
  (let [new-file (add-annotations file (first list) posneg)
        new-list (rest list)]
    (if (< 1 (count list))
      (recur new-file new-list posneg)
      file)))

(defn creat-one-gene [file list nontype]
  (let [new-file (add-annotations-gene file (first list) nontype)
        new-list (rest list)]
    (if (< 1 (count list))
      (recur new-file new-list nontype)
      file)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;OPPORATION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn gff3>dataset []
  (println "Opening a TSV file") 
  (def reference (ii/read-dataset "/home/yosh/datafiles/genes/merlin.gff3" :delim \tab))
  
  (println "Renameing colums")
  (def renamed (i/rename-cols
                {:col0     :seqid
                 :col1     :source
                 :col2     :type
                 :col3     :pt1
                 :col4     :pt2
                 :col5     :score
                 :col6     :strand
                 :col7     :phase
                 :col8     :attributes}
                reference))
  
  (def  attributed (add-id "gene" renamed))
  
  (def scrubed
    (->> attributed
         (i/$ [:seqid :gene :type :strand :pt1 :pt2 ])))
  (def empt       (empty-ref-db scrubed :type  235646))
  (def scrubed+   (i/$where {:strand {:$eq "+"}} scrubed))
  (def scrubed-   (i/$where {:strand {:$eq "-"}} scrubed))
  (def rescrubed  (i/$where {:type {:$eq "gene"}} scrubed))
  (def rescrubed+ (i/$where {:strand {:$eq "+"}} rescrubed))
  (def rescrubed- (i/$where {:strand {:$eq "-"}} rescrubed))
  (def ji+ (creat-one empt     (get-list scrubed+ :type) "+"))
  (def ji- (creat-one ji+      (get-list scrubed- :type) "-"))
  (def h+  (creat-one-gene ji- (get-list rescrubed+ :gene) :gene+))
  (def h+- (creat-one-gene h+  (get-list rescrubed- :gene) :gene-)))



    
