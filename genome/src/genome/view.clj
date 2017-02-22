(ns genome.view)
(require '[clojure.java.io :as io]
         '[incanter.core :as i]
         '[incanter.datasets :as id]
         '[incanter.io :as ii ]
         '[incanter.charts :as c]
         '[incanter.stats :as st]
         '[clojure.string :as s]
         '[clojure.data.csv :as csv]))

(def pied (ii/read-dataset file_out :header true))


(i/view (c/xy-plot
         :loc
         :pie
         :x-label "Position"
         :y-label "diversity"
         :title "02-519Pb AKA''Blip'"
                                        ;:legend true
         :data pied)))

(-> (c/xy-plot
     :loc
     :pie
     :x-label "Position"
     :y-label "diversity"
     :title "02-519Pb AKA''Blip'"
                                        ;:legend true
     :data pied) 
    (c/add-lines
     :loc
     (map #(float (/ % 100000)) (i/$ :cov pied)) ;:cov
     :data pied)
    (i/view))   

  (with-open [f-out (io/writer file_out)]
    (csv/write-csv f-out [(map name (i/col-names pied))])
    (csv/write-csv f-out (i/to-list pied))
