(ns genome.compare
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]
           [genome.database   :as gd ]
           [genome.stats      :as gs ]
           [genome.pop        :as p  ]
           [genome.consvar    :as cv ]
           [genome.dna2aa     :as da ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TESTING SNP TREND
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn snp-precent [snp file]
  "snp sould be the string of the keyword"
  (->> file
       (i/add-derived-column
        (keyword (str snp "fq"))
        [(keyword snp) :depth]
        #(if (= %2 0.0)
           0.0
           (/ %1 %2)))))



(defn add-snp-precent [file]
  (->> file
       (snp-precent "T")
       (snp-precent "C")
       (snp-precent "A")
       (snp-precent "G")))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE NEW COLUMN NAME
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn  col-rename [file ser_num]
  "This one is for adding a serial number at the end of a col name"
  (let [dont-change #{ :loc :merlin :gene+ :gene- :CDS+ :CDS- :exon+ :exon-}
        old-cols (apply vector (remove #(contains? dont-change %)
                                       (i/col-names file)))
        new-cols (->> old-cols
                      (map #(keyword (subs (str % ser_num) 1)))
                      (apply vector))]
    (->> new-cols
         (interleave old-cols)
         (apply assoc {}))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;UNITE TWO DATASET AT COMMON SITES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn add-row [file]
  "adds nil as first row because of the join bug"
  (let [ it_be (i/col-names file)]
    (i/conj-rows
     (i/dataset
      it_be
      [(vec (take (count it_be) (repeat nil)))])
     file)))

(defn merge-prep [file]
  "creates a dataset which contains all sites with allele frequency"
  (add-row (add-snp-precent file)))

(defn unite
  ([file1 file2]
   "Adds only file2 rows that have a common :loc value with file1"
   (let [set1   (merge-prep file1)
         set2   (merge-prep file2)
         coled1 (i/rename-cols (col-rename set1 1) set1)
         coled2 (i/rename-cols (col-rename set2 2) set2)]
     (->> coled1
          (i/$join [:loc :loc] coled2)
          (i/$where {:A2 {:$ne nil}}))))
  ([file1 file2 file3]
   (let [set1   (merge-prep file1)
         set2   (merge-prep file2)
         set3   (merge-prep file3)
         coled1 (i/rename-cols (col-rename set1 1) set1)
         coled2 (i/rename-cols (col-rename set2 2) set2)
         coled3 (i/rename-cols (col-rename set3 3) set3)]
     (->> coled1
          (i/$join [:loc :loc] coled2)
          (i/$join [:loc :loc] coled3)
          (i/$where {:A2 {:$ne nil}}))))
  ([file1 file2 file3 file4]
   (let [set1   (merge-prep file1)
         set2   (merge-prep file2)
         set3   (merge-prep file3)
         set4   (merge-prep file4)
         coled1 (i/rename-cols (col-rename set1 1) set1)
         coled2 (i/rename-cols (col-rename set2 2) set2)
         coled3 (i/rename-cols (col-rename set3 3) set3)
         coled4 (i/rename-cols (col-rename set4 4) set4)]
     (->> coled1
          (i/$join [:loc :loc] coled2)
          (i/$join [:loc :loc] coled3)
          (i/$join [:loc :loc] coled4)
          (i/$where {:A2 {:$ne nil}})))))

(def p-unite (memoize unite))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALCULATING MUTATION RATE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn mut-rate [file]
  "input.inc should be a synonimous mutation only"
  (let [com_syn+ (->> file
                      (i/$where {:CDS+ {:$ne "-" }})
                      (i/$where (i/$fn [pi1] (= pi1 0.0 )))
                      (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2] 
                                       (= majorf+1 minorf+1 majorf+2 minorf+2))))
        com_syn- (->> file
                      (i/$where {:CDS- {:$ne "-" }})
                      (i/$where (i/$fn [pi1] (= pi1 0.0 )))
                      (i/$where (i/$fn [majorf-1 minorf-1 majorf-2 minorf-2]
                                       (= majorf-1 minorf-1 majorf-2 minorf-2))))
        mut_syn+ (->> com_syn+
                      (i/$where (i/$fn [maj_p+1 maj_p+2 min_p+2]
                                       (or (not= maj_p+1 maj_p+2)
                                           (not= maj_p+1 min_p+2))))
                      i/nrow)
        mut_syn- (->> com_syn-
                     (i/$where (i/$fn [maj_p-1 maj_p-2 min_p-2]
                                      (or (not= maj_p-1 maj_p-2)
                                          (not= maj_p-1 min_p-2))))
                     i/nrow)
        com_sum  (double (+ (i/nrow com_syn+) (i/nrow com_syn-)))]

    (if (= 0.0 com_sum)
      0.0
      (/ (+ mut_syn+ mut_syn-) com_sum)))) 
    
       
  
    
