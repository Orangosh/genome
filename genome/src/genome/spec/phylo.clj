(ns genome.spec.phylo
  (require [clojure.java.io  :as io ]
           [clojure.string   :as str]
           [incanter.core    :as i  ]
           [incanter.io      :as ii ]
           [incanter.stats   :as st ]
           [incanter.charts  :as c  ]
           [clojure.data.csv :as csv]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; GET A TABLE WITH GENE START END FROM gff3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn pairs [type col]
  "parses through the ID line of gff3, creates a map returns a value"
  (->> col
       (re-seq  #"[^;]*")
       (map #(re-seq #"[^=]*" %))
       (map (fn [x] (filter #(not= "" %) x)))
       (filter #(= 2 (count %)))
       flatten
       (apply hash-map)
       (#(% type))))

(defn gene-list [file]
  "Takes gff3 file and derives gene list"
  (->> (ii/read-dataset file :delim \tab)
       (i/rename-cols {:col2 :type :col3 :pt1
                       :col4 :pt2  :col6 :strand
                       :col8 :attributes})
       (i/$where {:type {:eq "gene"}})
       (i/add-derived-column
        :gene
        [:attributes]
        #(pairs "gene" %))
       (i/$  [:gene :pt1 :pt2])
       (i/to-vect)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; PRINT THE GENE INTO FASTA FILE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-gene [[gene pt1 pt2] sample_name col file home] ;; add home
  (let [sequence (->> file
                      (i/$where (i/$fn [ref-loc]
                                       (and (<= pt1 ref-loc)
                                            (>= pt2 ref-loc))))
                      (i/$ col)
                      (remove #(= "-" %)) 
                      (apply str)
                      (partition-all 80)
                      (map #(apply str %)) )
        name     (str ">" sample_name)]
    (spit (str home "/phylo/input/" gene ".fas")
          (str name "\n" (apply str (map #(str % "\n") sequence)))
          :append true)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn fasta>inc [file]
  "file- the alinged ref + consensus sequences"
  (let [fasta (with-open [rdr (io/reader file)]
                (->>(line-seq rdr)
                    doall
                    (partition 2)))
        seq-col  (->> fasta
                      (map second)
                      (map #(vec (str/split % #"")))
                      (apply i/conj-cols)
                      (i/$where {:col-0 {:$ne "-"}}))
        name-col (->> fasta
                      (map #(keyword (subs (first %) 1)))
                       vec)
        locless  (i/rename-cols
                  (apply assoc {} (interleave (i/col-names seq-col) name-col))
                  seq-col)]
    (i/add-column :ref-loc (vec (range 1 (inc (i/nrow locless)))) locless)))


(defn inc>gene [gff-list file col home]
  (let [sample_name  (subs (str col) 1)]
    (map #(get-gene % sample_name col file home)
         gff-list)))


(defn gene>print [gff file home]
  "gff3- gff3 file"
  "file- The aligned ref + consensus sequences"
  "home- Target directory"
  (let [gff-list (gene-list gff )
        inc-file (fasta>inc file)]
    (map #(inc>gene gff-list inc-file % home)
         (drop-last (i/col-names inc-file)))))

