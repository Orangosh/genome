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

(def cov_un [:cov_un [:Tun :Aun :Cun :Gun]]);for variants calling
(def depth  [:depth  [:T   :A   :C   :G  ]]);after cleaning

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

(def maj_un [:maj_un+ [:Tun :Aun :Cun :Gun]]);for variants calling
(def major  [:maj+    [:T   :A   :C   :G  ]]);after minor allele variants

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
        (last  con_type) 
        #(get-major  %1 %2 %3 %4))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALL VARIANTS-MINOR ALLELE SEQUENCE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def min_un [:min_un+ [:Tun :Aun :Cun :Gun]]);for variants calling
(def minor  [:min+    [:T   :A   :C   :G  ]]);after minor allele variants

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
;;PRE-POISSON FILTERING ALLELE FREQUENCY PRECENTAGE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;                                                                                                  

(defn get-freq  [col_val maj_un cov_un]
  "Take a site read base and determent if it is in minimal minor allele def"
 (if (= "-"  maj_un)
   0.0
   (double (/ col_val cov_un))))


(defn add-freq-col [col_var col_name file]
  "Adds one column of nunleotide fithered by minor allele freq"
  (->> file
       (i/add-derived-column
        col_name
        [col_var :maj_un+ :cov_un]
	#(get-freq %1 %2 %3))))

(defn get-second [A T C G]
  "adds new base var column after iteration over one base column"
  (->> (vector A T C G)
       sort
       reverse
       second))


(defn minorify [file]
  "Adds four columns of nucleotide filthered by minor allele frequency"
  (->> file
       (add-freq-col :Aun :Atempfr)
       (add-freq-col :Tun :Ttempfr)
       (add-freq-col :Cun :Ctempfr)
       (add-freq-col :Gun :Gtempfr)
       (i/add-derived-column
        :tempminfr
        [:Atempfr :Ttempfr :Ctempfr :Gtempfr]
        #(get-second %1 %2 %3 %4))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ERROR FILTERING ASSUMING POISSON DISTRIBUTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn poisson-median [lambda col_val upper_case lower_case maj_un cov_un p]
  "Takes a site read base and using median calculates P asuming Poisson distribution"
 (if (= "-"  maj_un)
   0
   (let [cov_lambda (* cov_un lambda)]
     (if (and (or (<= 1 cov_lambda);if the coverage corrected lambda is < = than it means that the coverage is too low
                  (not= 1 col_val))
              (< p (st/cdf-poisson col_val :lambda cov_lambda))
              (and (> upper_case 1) (> lower_case 1))) ;make sure thate the varient exists at both strands
       col_val
       0))))


;adds new base var column after iteration over one base column
(defn pois-correct-median [lambda col_var upper_case lower_case col_name p file]
  "Adds one column of nunleotide fithered by poison distribution"
  (->> file
       (i/add-derived-column
        col_name
        [col_var upper_case lower_case :maj_un+ :cov_un]
        #(poisson-median lambda %1 %2 %3 %4 %5 p))))

(defn minorify-post-poisson [file]
  "Adds four columns of nucleotide filthered by minor allele frequency after poison"
  (->> file
       (add-freq-col :A :Afr)
       (add-freq-col :T :Tfr)
       (add-freq-col :C :Cfr)
       (add-freq-col :G :Gfr)
       (i/add-derived-column
        :minfr
        [:Afr :Tfr :Cfr :Gfr]
        #(get-second %1 %2 %3 %4))))

(defn poissonize-median [p_value file]
  "Adds four columns of nucleotide filthered by poison distribution"
  (let [p      (- 1 p_value)
        lambda (st/median (i/$ :tempminfr (i/$where (i/$fn [tempminfr] (> tempminfr 0.0)) file)))]
    (->> file
         (pois-correct-median lambda :Aun \A \a :A p)
         (pois-correct-median lambda :Tun \T \t :T p)
         (pois-correct-median lambda :Cun \C \c :C p)
         (pois-correct-median lambda :Gun \G \g :G p)
         (minorify-post-poisson))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;ERROR FILTERING METHOD
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-minor-allele [pois_p file]
  "first value the p-value of poisson if 1- no possion filtering"
  (if (= 1 pois_p)
    (i/rename-cols {:Aun     :A   :Tun     :T   :Cun     :C   :Gun     :G
                    :Atempfr :Afr :Ttempfr :Tfr :Ctempfr :Cfr :Gtempfr :Gfr :tempminfr :minfr} (minorify file))
    (poissonize-median pois_p (minorify file))))
 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                         !!!!! DEPRACETED !!!!!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ERROR FILTERING ASSUMING POISSON DISTRIBUTION- NOT USED
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(def T_matrix
  "Transition matrix as per illumina reads for varent calling."
  (i/dataset [:major  :Aun  :Tun :Cun  :Gun] 
             [["A"        0   0.92  1.75  2.00] 
              ["T"      0.94   0    2.41  1.82] 
              ["C"      1.81  1.19    0   1.78]
              ["G"      1.21  1.15  1.81    0 ]]))


(defn poisson-T [col_var col_val maj_un cov_un p]
  "Takes a site read base and Using T_matrix calculates P asuming Poisson distribution"
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
(defn pois-correct-T [col_var col_name p file]
  "Adds one column of nunleotide fithered by poison distribution"
  (->> file
       (i/add-derived-column
        col_name
        [col_var :maj_un+ :cov_un]
        #(poisson-T col_var %1 %2 %3 p))))

(defn poissonize-T [p_value file]
  "Adds four columns of nucleotide filthered by poison distribution"
  (let [p (- 1 p_value)]
    (->> file
         (pois-correct-T :Aun :A p)
         (pois-correct-T :Tun :T p)
         (pois-correct-T :Cun :C p)
         (pois-correct-T :Gun :G p))))












