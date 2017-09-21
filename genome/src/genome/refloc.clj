(ns genome.refloc
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [genome.database   :as gd ]
           [clojure.data.csv  :as csv]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; This one takes only the consensus sequence from the multiple alingment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

(defn get-refinc[fasta_file]
  "/home/yosh/datafiles/Consensus/CMVconsensus/exphcmv.fasta"
  (let [ fasta (with-open [rdr (clojure.java.io/reader fasta_file)]
                 (reduce conj [] (line-seq rdr)))
        ref-name (first fasta)
        ref-seq  (second fasta)
        ref-loc  (get-loc-vec ref-seq [] 1)
        ]
    (->> (i/conj-cols (vec ref-seq) (vec (range 1 (+ 1 (count ref-seq)))) ref-loc)
         (i/rename-cols {:col-0 :ref-seq :col-1 :loc :col-2 :ref-loc} ))))

(def m-get-refinc (memoize get-refinc))

(defn refer-ann [refset file]
  "/mnt/data/ebv/consensus/ebv-ref-NC_007605.1.fasta"
  (let [ann_ref (m-get-refinc refset)]
    (i/$join [:loc :loc] ann_ref file)))
 
