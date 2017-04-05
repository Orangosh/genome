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
;TESTING SNP TREND
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn calc-coved [file]
  "Calculation corrected coverage"
  (->> file
       (i/add-derived-column
        :p_cov
        [:Tpois :Apois :Gpois :Cpois]
        #(+ %1 %2 %3 %4))))

(defn snp-precent [snp file]
  "snp sould be the string of the keyword"
  (->> file
       (i/add-derived-column
        (keyword (str snp "-per"))
        [(keyword snp) :p_cov]
        #(if (= %2 0)
           0.0
           (double (/ %1 %2))))))

(defn add-snp-precent [file]
  (->> file
       (snp-precent "Tpois")
       (snp-precent "Cpois")
       (snp-precent "Apois")
       (snp-precent "Gpois")))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE NEW COLUMN NAME
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def sample1 {:Tpois :T1 :Tpois-per :Tper1
              :Cpois :C1 :Cpois-per :Cper1
              :Apois :A1 :Apois-per :Aper1
              :Gpois :G1 :Gpois-per :Gper1
              :p_cov :cov1 :pi_pois :pi1})

(def sample2 {:Tpois :T2 :Tpois-per :Tper2
              :Cpois :C2 :Cpois-per :Cper2
              :Apois :A2 :Apois-per :Aper2
              :Gpois :G2 :Gpois-per :Gper2
              :p_cov :cov2 :pi_pois :pi2})


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;UNITE TWO DATASET AT COMMON SITES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn unite [file1 file2]
  "Adds only file2 rows that have a common :loc value with file1"
  (let [coled1 (i/rename-cols sample1 file1)
        coled2 (i/rename-cols sample2 file2)]
    (i/$where {:A2 {:$ne nil}}
              (i/$join [:loc :loc] coled2 coled1))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn create [file1 file2]
  (let [p_covd1 (calc-coved file1)
        p_covd2 (calc-coved file2)
        snp1   (add-snp-precent p_covd1)
        snp2   (add-snp-precent p_covd1)]
    (->>(unite snp1 snp1)
        (i/$ [:loc :cov1
              :A1 :Aper1 :T1 :Tper1
              :G1 :Gper1 :C1 :Cper1 :pi1
              :cov2
              :A2 :Aper2 :T2 :Tper2
              :G2 :Gper2 :C2 :Cper2 :pi2])))) 

(def comon_pied (i/$where {:pi1 {:$ne 0.0}} (i/$where {:pi2 {:$ne 0.0}} b)))

;(def c (i/$ [:cov1 :Tpois :Tpois-precent  :Cpois :Cpois-precent  :Apois :Apois-precent  :Gpois :Gpois-precent] b))
  
