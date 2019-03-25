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
  (def mean_cov (st/mean (i/$ :depth file)))
  (def nucleotide_diversity (/  (i/sum (i/$ :pi file)) (i/nrow file)))
  (def segregating_sites (count (filter #(< 0 %) (i/$ :pi file))))
  (def binned (p/bin-sfs 10 file))
  (println "mean coverage " mean_cov)
  (println "Nucleotide diversity: " nucleotide_diversity)
  (println "Segregating Sites: " segregating_sites)
  (println "folded and binned (10) SFS " binned))


;; Make sure incstates.clj fn are loaded!!!

 
(defn get-all-mean [samples names]
  (st/mean 
       (map (fn [y] (st/mean (i/$ :pi y)))
            (map (fn [z] (z samples))
                 (map (fn [w] (keyword w))
                      names)))))



(get-all-mean ebv_samples ebv_primaries_names)
(get-all-mean hhv6_samples hhv6_primaries_names)
(get-all-mean hcmv_samples hcmv_primaries_names)

