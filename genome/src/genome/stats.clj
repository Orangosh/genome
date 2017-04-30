(ns genome.stats
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.datasets :as id]
           [incanter.io :as ii ]
           [incanter.charts :as c]
           [incanter.stats :as st]
           [clojure.string :as s]
           [clojure.data.csv :as csv]
           [incanter.distributions :as dst]
           [incanter.zoo :as z]
           [genome.pop :as p]))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CLEAN ROWS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn row-clean [file column to_clean]
  (->> file
       (i/$where {column {:$ne to_clean} })))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SUMMARY STATISTICS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn stat-report [file]
  (def mean_cov (st/mean (i/$ :cov_p file)))
  (def nucleotide_diversity (/  (i/sum (i/$ :pi_pois file)) (i/nrow file)))
  (def segregating_sites (count (filter #(< 0 %) (i/$ :pi_pois file))))
  (def binned (p/bin 10 file))
  (println "mean coverage " mean_cov)
  (println "Nucleotide diversity: " nucleotide_diversity)
  (println "Segregating Sites: " segregating_sites)
  (println "folded and binned (10) SFS " binned))


;(def incar (ii/read-dataset "/home/yosh/datafiles/incanted" :header true))

