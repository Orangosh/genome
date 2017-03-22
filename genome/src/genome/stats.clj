(ns genome.stats
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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DIVERSITY
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Calculats pi for each row
(defn pi [cov T A C G] 
  (if (>=  cov 2) 
    (double (/(+ (* T A) (* T C)
                 (* T G) (* A C) 
                 (* A G) (* G C))
              (/ (* cov (- cov 1))
                 2)))
    0))

(def pi_un [:pi_un [:c_cov :Tun :Aun :Cun :Gun]]);for variants calling
(def pi_pois [:pi_pois [:c_cov :Tpois :Apois :Cpois :Gpois]]);after variants

(defn pise [file pi_type]
  (->> file
       (i/add-derived-column
         (first pi_type)
         (last pi_type) 
         #(pi %1 %2 %3 %4 %5))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SLIDING WINDOW
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Creates a sliding window from a column :pi in file and adds a new col :pislide 
(defn slide [file scanned_column win_size]
  (i/add-column
   :sliding
   (->> (i/$ scanned_column file)
        (partition win_size 1)
        (map #(/ (apply + %) win_size))
        (concat (take (dec win_size) (repeat 0))))
   file))

;Very similar with a built in function
(defn glide [file scanned_column win_size]
  (i/add-column
   :gliding
   (->> (i/$ scanned_column file)
        (z/roll-mean win_size) ;instead of roll-min we can use roll-apply func
        (concat (take (dec win_size) (repeat 0))))
   file))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALL CONCENSUS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def consus_un [:consus_un [:Tun :Aun :Cun :Gun]]);for variants calling
(def consus_pois [:consus_pois [:Tpois :Apois :Cpois :Gpois]]);after variants

(defn get-max [T A C G] 
  (let [get_map {"T" T "A" A "C" C "G" G}] 
    (->> get_map 
         (keep #(when (= (val %) (apply max (vals get_map))) (key %))) 
         rand-nth))) ;rand-nth will choose equally apearing nucleotide at a site

(defn concensus [file con_type]
  (->> file
       (i/add-derived-column
        (first con_type)
        (last con_type) 
        #(get-max  %1 %2 %3 %4))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALL VARIANTS ASSUMING POISSON DISTRIBUTION
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
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ALLELE FREQUENCY SPECTRA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn unfolded-SFS [ref T A C G] 
  (if-not (= ref "-")
    (let [f { "T" T "A" A "C" C "G" G}] 
      (apply max (filter #(not= (f ref) %) [ T A C G])))
    "-")) 

;Calculates unfolded site allele frequency for all minor alleles
(defn multi-SFS [ref T A C G] 
  (if-not (= ref "-")
    (let [f { "T" T "A" A "C" C "G" G}] 
      (i/sum (filter #(not= (f ref) %) [ T A C G])))
    "-"))

(defn folded-SFS [mean_cov ref c_cov T A C G] 
  (let [f { "T" T "A" A "C" C "G" G}] 
      (second (reverse (sort [ T A C G])))))

(defn folded-SFS [mean_cov ref c_cov T A C G] 
  (->> [T A C G]
       sort
       reverse
       second
       (* (/ mean_cov c_cov))
       double))

(defn SFS [file SFS-type]
  (let [mean_cov (st/mean (i/$ :c_cov file))]
    (->> file
         (i/add-derived-column
          :sfs
          [:ref :c_cov :Tpois :Apois :Cpois :Gpois]
          #(SFS-type mean_cov %1 %2 %3 %4 %5 %6)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;BINNING
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn bin [n-bins file] ;bin range 0-14 into 5 bins (bin 5 (range 15))
  (let [bin-array (i/$ :sfs file)
        min-freq (apply min bin-array)
        max-freq (apply max bin-array)
        range-freq (- max-freq min-freq)
        bin-fn (fn [freq](-> freq
                          (- min-freq)
                          (/ range-freq)
                          (* n-bins)
                          (int)
                          (min (dec n-bins))))]
    (->> (map bin-fn bin-array)
         frequencies
         sort)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CLEAN ROWS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn row-clean [file column to_clean]
  (->> file
       (i/$where {column {:$ne to_clean} })))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SUMMARY STATISTICS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn stat-report [file]
  (def nucleotide_diversity (/  (i/sum (i/$ :pi_pois file)) (i/nrow file)))
  (def segregating_sites (count (filter #(< 0 %) (i/$ :pi_pois file))))
  (def binned (bin 10 file))
  (println "Nucleotide diversity: " nucleotide_diversity)
  (println "Segregating Sites: " segregating_sites)
  (println "folded and binned (10) SFS " binned))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TESTING PIPELINE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



(defn statistics [file]
;"/home/yosh/datafiles/incar"
;(def incar (ii/read-dataset file :header true))
  (println "Calculation corrected coverage")
  (def c_cov (->> file
                  (i/add-derived-column
                   :c_cov
                   [:Tun :Aun :Gun :Cun]
                   #(+ %1 %2 %3 %4))))
  (println "Creating first concensus")
  (def conded (concensus c_cov consus_un))
  (println "Correcting read errors")
  (def pois (poissonize 0.05 conded))
  (def scrubed (i/$ [:r_seq :loc :ref :consus_un :cov :c_cov
                     :Tpois :Apois :Cpois :Gpois] pois))
  (println "Calculating nucleotide diversity")
  (def pied (pise scrubed pi_pois))
  (println "Calculating folded allele frequency spectra")
  (def sfsd (SFS pied folded-SFS))
  (def row_cleaned (row-clean sfsd :ref "-"))
  (println "SUMMARY STATISTICS:")
  (stat-report sfsd))


  
