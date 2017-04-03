(ns genome.compare
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.datasets :as id]
           [incanter.io :as ii ]
           [incanter.charts :as c]
           [incanter.stats :as st]
           [clojure.string :as s]
           [clojure.data.csv :as csv]
           [genome.database :as gd]
           [genome.stats :as gs]
           [genome.pop :as p]
           [genome.consvar :as cv]
           [genome.dna2aa :as da]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;UNITING 2 FILES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(def sample1  {:Tpois :T1
               :Cpois :C1
               :Apois :A1
               :Gpois :G1
               :pi_pois :pi1})

(def sample2  {:Tpois :T2
               :Cpois :C2
               :Apois :A2
               :Gpois :G2
               :pi_pois :pi2})

(def sample3  {:Tpois :T3
               :Cpois :C3
               :Apois :A3
               :Gpois :G3
               :pi_pois :pi3})

(defn unite [file1 file2]
  "Adds only file2 rows that have a common :loc value with file1"
  (let [coled1 (i/rename-cols sample1 file1)
        coled2 (i/rename-cols sample2 file2)]
    (i/$where {:A2 {:$ne nil}}
              (i/$join [:loc :loc] coled2 coled1))))



(def b (i/$ [:loc :pi1 :pi2
             :A1 :T1 :G1 :C1
             :A2 :T2 :G2 :C2] a)) 
(def nu (i/$where {:pi1 {:$ne 0.0}} (i/$where {:pi2 {:$ne 0.0}} b)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TESTING PIPELINE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
