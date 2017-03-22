(ns genome.core
  (:gen-class)
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.datasets :as id]
           [incanter.io :as ii ]
           [incanter.charts :as c]
           [incanter.stats :as st]
           [clojure.string :as s]
           [clojure.data.csv :as csv]
           [genome.database :as gd]
           [genome.stats :as gs]))

(defn -main [file_in file_out]
  (println "Wellcome to clojure- starting incanter")
  (gd/create-db file_in)
  (gs/statistics gd/finalized)
  (with-open [f-out (io/writer file_out)]
    (csv/write-csv f-out [(map name (i/col-names gs/row_cleaned))])
    (csv/write-csv f-out (i/to-list gs/row_cleaned))))


 
