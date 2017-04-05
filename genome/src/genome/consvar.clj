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
;CALL CONSENSUS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def consus_un [:consus_un [:Tun :Aun :Cun :Gun]]);for variants calling
(def consus_pois [:consus_pois [:Tpois :Apois :Cpois :Gpois]]);after variants

(defn get-max [T A C G] 
  (let [get_map {"T" T "A" A "C" C "G" G}] 
    (->> get_map 
         (keep #(when (= (val %) (apply max (vals get_map))) (key %))) 
         rand-nth))) ;rand-nth will choose equally apearing nucleotide at a site

(defn consensus [file con_type]
  (->> file
       (i/add-derived-column
        (first con_type)
        (last con_type) 
        #(get-max  %1 %2 %3 %4))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ERROR FILTERING ASSUMING POISSON DISTRIBUTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;transition matrix as per illumina reads for varent calling.
(def T_matrix (i/dataset [:consus_un    :Aun  :Tun :Cun  :Gun] 
                         [["A"         0   0.92 1.75  2.00] 
                          ["T"       0.94   0   2.41  1.82] 
                          ["C"       1.81  1.19   0   1.78]
                          ["G"       1.21  1.15 1.81    0 ]]))

;take a site read base and calculat P asuming Poisson distribution
(defn poisson [col_var col_val consus c_cov p]
  (let [lambda (->> T_matrix
                    (i/$where {:consus_un {:$eq consus}})
                    (i/$ col_var)
                    (* (/ c_cov 1000)))]
    (if (< p (st/cdf-poisson col_val :lambda lambda))
      col_val
      0)))
  
;adds new base var column after iteration over one base column
(defn pois-correct [col_var col_name p file]
  (->> file
       (i/add-derived-column
        col_name
        [col_var :consus_un :c_cov]
        #(poisson col_var %1 %2 %3 p))))

;add four columns
(defn poissonize [p_value file] 
  (let [p (- 1 p_value)]
    (->> (pois-correct :Aun :Apois p file)
         (pois-correct :Tun :Tpois p)
         (pois-correct :Cun :Cpois p)
         (pois-correct :Gun :Gpois p))))

(defn calc-coved [file]
  "Calculation corrected coverage"
  (->> file
       (i/add-derived-column
        :p_cov
        [:Tun :Aun :Gun :Cun]
        #(+ %1 %2 %3 %4))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALL VARIANTS 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def var_un [:var_un [:Tun :Aun :Cun :Gun]]);for variants calling
(def var_pois [:var_pois [:Tpois :Apois :Cpois :Gpois]]);after variants

;gets the variant amino acid

(defn get-var-nuc [T A C G] 
                (let [get_map {"T" T "A" A "C" C "G" G}
                      minor_allele (second (reverse (sort (vals get_map))))]
                  (cond
                    (= minor_allele 0) "-"
                    :else (->> get_map 
                               (keep #(when (= (val %) minor_allele) (key %)))
;rand-nth will choose equally apearing nucleotide at a site
                               rand-nth))))                

(defn minor-alleles [file con_type]
  (->> file
       (i/add-derived-column
        (first con_type)
        (last con_type) 
        #(get-var-nuc  %1 %2 %3 %4))))
