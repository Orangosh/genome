(ns genome.consvar
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.datasets :as id]
           [incanter.io :as ii ]
           [incanter.charts :as c]
           [incanter.stats :as st]
           [clojure.string :as s]
           [clojure.data.csv :as csv]
           [incanter.distributions :as dst]
           [incanter.zoo :as z]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CORRECT COVERAGE 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def cov_un [:cov_un [:Tun   :Aun   :Cun   :Gun  ]]);for variants calling
(def cov_p  [:cov_p  [:Tpois :Apois :Cpois :Gpois]]);after poison variants
(def depth  [:depth  [:T     :A     :C     :G  ]]);after cleaning

(defn calc-coved [con_type file]
  "Calculation corrected coverage"
  (->> file
       (i/add-derived-column
        (first con_type)
        (last con_type) 
        #(double (+ %1 %2 %3 %4)))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALL CONSENSUS- MAJOR ALLELE SEQUENCE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def maj_un [:maj_un+ [:Tun   :Aun   :Cun   :Gun  ]]);for variants calling
(def maj_p  [:maj_p+  [:Tpois :Apois :Cpois :Gpois]]);after variants
(def maj  [:maj+  [:T     :A     :C     :G ]]);after minor allele variants

(defn get-major [T A C G] 
  (let [get_map {"T" T "A" A "C" C "G" G}] 
    (if (= 0.0 (double (+ T A C G)))
      "-"
      (->> get_map 
           (keep #(when (= (val %) (apply max (vals get_map))) (key %)))
           rand-nth)))) ;rand-nth will choose equally apearing nucleotide at a site

(defn major-allele [con_type file]
  (->> file
       (i/add-derived-column
        (first con_type)
        (last con_type) 
        #(get-major  %1 %2 %3 %4))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALL VARIANTS-MINOR ALLELE SEQUENCE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def min_un [:min_un+ [:Tun   :Aun   :Cun   :Gun  ]]) ;for variants calling
(def min_p  [:min_p+  [:Tpois :Apois :Cpois :Gpois]]);after variants
(def min    [:min+    [:T  :A  :C  :G ]]);after minor allele variants

(defn get-minor [T A C G] 
  "Gets minor allels rand-nth- chooses equally apearing nucleotide at a site"
  (let [get_map   {"T" T "A" A "C" C "G" G}
        minor_allele (second (reverse (sort (vals get_map))))]
    (if (= 0.0 (double (+ T A C G)))
      "-"
      (cond
        (= minor_allele 0)
        (->> get_map 
             (keep #(when (= (val %) (apply max (vals get_map))) (key %))) 
             rand-nth)
        :else
        (->> get_map 
             (keep #(when (= (val %) minor_allele) (key %)))
             rand-nth)))))           

(defn minor-allele [con_type file]
  (->> file
       (i/add-derived-column
        (first con_type)
        (last con_type) 
        #(get-minor  %1 %2 %3 %4))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ERROR FILTERING ASSUMING POISSON DISTRIBUTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(def T_matrix
  "Transition matrix as per illumina reads for varent calling."
  (i/dataset [:major  :Aun  :Tun :Cun  :Gun] 
             [["A"        0   0.92  1.75  2.00] 
              ["T"      0.94   0    2.41  1.82] 
              ["C"      1.81  1.19    0   1.78]
              ["G"      1.21  1.15  1.81    0 ]]))


(defn poisson [col_var col_val maj_un cov_un p]
  "Take a site read base and calculat P asuming Poisson distribution"
 (if (= "-"  maj_un)
   0
   (let [lambda (->> T_matrix
                     (i/$where {:major {:$eq maj_un}})
                     (i/$ col_var)
                     (* (/ cov_un 1000)))]
     (if (< p (st/cdf-poisson col_val :lambda lambda))
       col_val
       0))))
  
;adds new base var column after iteration over one base column
(defn pois-correct [col_var col_name p file]
  "Adds one column of nunleotide fithered by poison distribution"
  (->> file
       (i/add-derived-column
        col_name
        [col_var :maj_un+ :cov_un]
        #(poisson col_var %1 %2 %3 p))))


(defn poissonize [p_value file]
  "Adds four columns of nucleotide filthered by poison distribution"
  (let [p (- 1 p_value)]
    (->> file
         (pois-correct :Aun :A p)
         (pois-correct :Tun :T p)
         (pois-correct :Cun :C p)
         (pois-correct :Gun :G p))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ERROR FILTERING ASSUMING ALLELE FREQUENCY PRECENTAGE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn allele?  [col_var col_val maj_un cov_un minor_freq]
  "Take a site read base and determent if it is in minimal minor allele def"
 (if (= "-"  maj_un)
   0
   (if (<  minor_freq
           (double (/ col_val cov_un)))
     col_val
     0)))
  
;adds new base var column after iteration over one base column
(defn minor-correct [col_var col_name minor_freq file]
  "Adds one column of nunleotide fithered by minor allele freq"
  (->> file
       (i/add-derived-column
        col_name
        [col_var :maj_un+ :cov_un]
        #(allele? col_var %1 %2 %3 minor_freq))))


(defn minorallize [minor_freq file]
  "Adds four columns of nucleotide filthered by minor allele frequency"
  (->> file
       (minor-correct :Aun :A minor_freq)
       (minor-correct :Tun :T minor_freq)
       (minor-correct :Cun :C minor_freq)
       (minor-correct :Gun :G minor_freq)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ERROR FILTERING METHOD
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-minor-allele [method value file]
  (case method
    "minor allele" (minorallize value file)
    "poisson dist" (poissonize  value file)))
