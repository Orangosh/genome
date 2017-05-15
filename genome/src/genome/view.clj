(ns genome.view 
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]))

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
         :data file)))
  ([file cov_div]
  (-> (c/xy-plot
       :loc
       :pie
       :x-label "Position"
       :y-label "diversity"
       :title "02-519Pb AKA''Blip'"
      ;:legend true
       :data file) 
      (c/add-lines
       :loc
       (map #(float (/ % cov_div)) (i/$ :cov file)) ;:cov
       :data file)
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
        (i/view)))

  ([column file1 file2 file3 file4 file5 file6]
    (-> (c/xy-plot
         :loc
         column
         :x-label "frequency"
         :y-label "Location"
         :title   "Change in nucelotide diversity per site between 2 samples"
         :data file1) 
        (c/add-lines
         :loc
         column
         :data file2)
        (c/add-lines
         :loc
         column
         :data file3)    
        (c/add-lines
         :loc
         column
         :data file4)    
        (c/add-lines
         :loc
         column
         :data file5)    
        (c/add-lines
         :loc
         column
         :data file6)    
        (i/view))))

(defn view-gene [file col]
  (i/with-data (->> file
                    (i/$ col)
                    frequencies
                    (map vec)
                    vec
                    (i/dataset[col :sum]))
    (i/view (c/bar-chart col :sum))))

(defn view-gene [file col]
  (-> (c/bar-chart
       col
       :sum
       :data (->> file
                  (i/$ col)
                  frequencies
                  (map vec)
                  vec
                  (i/dataset[col :sum])))
      i/view))

(defn view-gene [file col file2 col2]
  (-> (c/bar-chart
       col
       :sum
       :data (->> file
                  (i/$ col)
                  frequencies
                  (map vec)
                  vec
                  (i/dataset[col :sum])))
      (c/add-lines
       col2
       :sum
       :data (->> file2
                  (i/$ col)
                  frequencies
                  (map vec)
                  vec
                  (i/dataset[col :sum])))
      i/view))
