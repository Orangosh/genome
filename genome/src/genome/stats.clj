(ns genome.stats
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.datasets :as id]
           [incanter.io :as ii ]
           [incanter.charts :as c]
           [incanter.stats :as st]
           [clojure.string :as s]
           [clojure.data.csv :as csv]))


(defn stat-report [file_in]
  (def pied (ii/read-dataset file_in :header true))
  (def nucleotide_diversity (/  (i/sum (i/$ :pie pied)) (i/nrow pied)))
  (def segregation_sites (count (filter #(< 0 %) (i/$ :pie pied))))
  (println "Nucleotide diversity: " nucleotide_diversity)
  (println "Segregating Sites: " segregation_site))

;creates a sliding window from a column :pi in file and adds a new col :pislide 
(defn win-slide [file win_size]
  (let [winset (ii/read-dataset file :header true)]
    (i/add-column
     :pislide
     (->> (i/$ :pi winset)
          (partition win_size 1)
          (map #(/ (apply + %) win_size))
          (concat (take (dec win_size) (repeat 0))))
     winset)))


(defn pi [T A G C] 
  (let [cov (+ T A G C)]
    (if (>=  cov 2) 
      (double (/(+ (* T A) (* T G)
                   (* T C) (* A G) 
                   (* A C) (* G C))
                (/ (* cov (- cov 1))
                   2)))
      0)))


(defn SFS [ref T A C G]
  (let [f { "T" T "A" A "C" C "G" G}] 
    (i/sum (filter #(not= (f ref) %) [ T A C G]))))


