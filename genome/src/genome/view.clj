(ns genome.view)
(require '[clojure.java.io :as io]
         '[incanter.core :as i]
         '[incanter.datasets :as id]
         '[incanter.io :as ii ]
         '[incanter.charts :as c]
         '[incanter.stats :as st]
         '[clojure.string :as s]
         '[clojure.data.csv :as csv]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;VIEW ONE FILE PROPERTIES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn pi-view
  ([file]
  (i/view (c/xy-plot
         :loc
         :pie
         :x-label "Position"
         :y-label "diversity"
         :title "Pi win size"
        ;:legend true
         :data pied)))
  ([file cov_div]
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
       (map #(float (/ % cov_div)) (i/$ :cov pied)) ;:cov
       :data pied)
      (i/view))))   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;WIEW COMPARRISSONES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn look
  "create a graph of the value comparison"
  ([column file]
   (i/view (c/xy-plot
            :loc
            column
            :x-label "frequency"
            :y-label "Location"
            :title   "Change in nucelotide diversity per site between 2 samples"
;:legend true
            :data file)))

  ([column early>inter early>late]
   (-> (c/xy-plot
        :loc
        column
        :x-label "frequency"
        :y-label "Location"
        :title   "Change in nucelotide diversity per site between 2 samples"
        :data early>inter) 
       (c/add-lines
        :loc
        column
        :data early>late)
       (i/view)))

  ([column early>inter early>late early>mom_sib]
    (-> (c/xy-plot
         :loc
         column
         :x-label "frequency"
         :y-label "Location"
         :title   "Change in nucelotide diversity per site between 2 samples"
         :data early>inter) 
        (c/add-lines
         :loc
         column
         :data early>late)
        (c/add-lines
         :loc
         column
         :data early>mom_sib)    
        (i/view))))
