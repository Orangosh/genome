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



(defn win-slide [file win_size]
  (concat (take (dec win_size) (repeat 0)) 
          (map #(/ (apply + %) win_size) 
               (partition win_size 1 file)))
  )

(def r (i/add-column :win_s (win-slide q 100) p))

