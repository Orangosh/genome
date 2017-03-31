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
         :title "Pi win size"
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
  (csv/write-csv f-out (i/to-list pied)))

(take-nth 2 (rest ( flatten (sort (frequencies (i/$ :sfs SFSd)))))) ;get sfs

(i/dataset [ :sfs :count] ;get a dataset wit sfs (from databas SFSd
           (vec (sort (frequencies (i/$ :sfs SFSd)))))

(i/view (c/xy-plot ;view sfs
         :sfs
         :count
         :x-label "Position"
         :y-label "diversity"
         :title "02-519Pb AKA''Blip'"
                                        ;:legend true
         :data (i/$ (range 0 20) :all j)))


(def idf2 (i/to-dataset (vec gs/binned2)))
(def idf (i/to-dataset (vec gs/binned)))

(i/view (c/xy-plot
         :col-0
         :col-1
         :x-label "frequency"
         :y-label "occurance"
         :title "SFS bin 100"
                                        ;:legend true
         :data idf2))

(i/view (c/xy-plot
         :col-0
         :col-1
         :x-label "frequency"
         :y-label "occurance"
         :title "SFS bin 10"
                                        ;:legend true
         :data idf))
