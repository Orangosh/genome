(ns genome.conloc
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]))


(defn get-conset[fasta_file]
  "/home/yosh/datafiles/Consensus/CMVconsensus/exphcmv.fasta"
  (let [fasta (with-open [rdr (clojure.java.io/reader fasta_file)]
                (reduce conj [] (line-seq rdr)))]
    (apply i/conj-cols (take-nth 2 (rest (map vec fasta))))))

(defn get-loc-vec [ref loc-vec loc]
  (if (< 0 (count ref))
    (if (= (first ref) \-)
      (recur
       (rest ref)
       (conj loc-vec loc)
       (loc))
      (recur
       (rest ref)
       (conj loc-vec loc)
       (inc loc)))
    loc-vec))

