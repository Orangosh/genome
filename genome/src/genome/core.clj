(ns genome.core
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.datasets :as id]
           [incanter.io :as ii ]
           [incanter.charts :as c]
           [incanter.stats :as st]
           [clojure.string :as s]
           [clojure.data.csv :as csv]
           [genome.database :as d]))

(defn -main   [file_in file_out]
  (println "Good morning Vietnam")
  (d/create-db file_in file_out)
  
  (println "opening")
  )
  


 
