(ns genome.conloc
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [genome.database   :as gd ]
           [clojure.data.csv  :as csv]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TESTING PIPELINE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
(defn get-conset[fasta_file]
  "/home/yosh/datafiles/Consensus/CMVconsensus/exphcmv.fasta"
  (let [fasta (with-open [rdr (clojure.java.io/reader fasta_file)]
                (reduce conj [] (line-seq rdr)))
        temp  (apply i/conj-cols (take-nth 2 (rest (map vec fasta))))]
    (->> (i/rename-cols (zipmap (i/col-names temp)
                                (vec (map keyword
                                          (map #(subs % 4)
                                               (take-nth 2 fasta)))))
                        temp)
         (i/rename-cols {:lin :merlin}))))

(defn get-loc-vec [ref loc-vec loc]
  (if (< 0 (count ref))
    (if (= (first ref) \-)
      (recur (rest ref)
             (conj loc-vec loc)
             loc)
      (recur (rest ref)
             (conj loc-vec loc)
             (inc loc)))
    loc-vec))

(def m-get-loc-vec (memoize get-loc-vec))

(defn add-count-cols [fasta_file]
  (let [file    (get-conset fasta_file)
        ref     (i/$ :merlin file)
        ref_loc (get-loc-vec ref [] 1)]
    (->> (i/add-column :ref-loc ref_loc file)
         (i/add-column :loc (range (i/nrow file))))))

(def m-add-count-cols (memoize add-count-cols))

(def refset (m-add-count-cols "/home/yosh/datafiles/Consensus/CMVconsensus/exphcmv.fasta"))

(def refset.inc "/home/yosh/datafiles/Consensus/CMVconsensus/refset.inc")

(with-open [f-out (io/writer refset.inc)]
    (csv/write-csv f-out [(map name (i/col-names refset))])
    (csv/write-csv f-out (i/to-list refset)))
