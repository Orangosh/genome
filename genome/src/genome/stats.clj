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


;in the future I would like to move all stats here from database

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DIVERSITY
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Calculats pi for each row
(defn pi [T A G C] 
  (let [cov (+ T A G C)]
    (if (>=  cov 2) 
      (double (/(+ (* T A) (* T G)
                   (* T C) (* A G) 
                   (* A C) (* G C))
                (/ (* cov (- cov 1))
                   2)))
      0)))

(defn pied [file]
  (->> file
       (i/add-derived-column
        :pi
        [:Tun :Aun :Gun :Cun]
        #(pi %1 %2 %3 %4))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SLIDING WINDOW
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Creates a sliding window from a column :pi in file and adds a new col :pislide 
(defn slide [file  scanned_column win_size]
  (let [winset (ii/read-dataset file :header true)]
    (i/add-column
     :sliding
     (->> (i/$ scanned_column winset)
          (partition win_size 1)
          (map #(/ (apply + %) win_size))
          (concat (take (dec win_size) (repeat 0))))
     winset)))

;Very similar with a built in function
(defn slide-try [file scanned_column win_size]
  (let [winset (ii/read-dataset file :header true)]
    (i/add-column
     :slidingTry
     (->> (i/$ scanned_column winset)
          (z/roll-mean win_size) ;instead of roll-min we can use roll-apply func
          (concat (take (dec win_size) (repeat 0))))
     winset)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALL CONCENSUS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-max [T A C G] 
  (let [get_map {"T" T "A" A "C" C "G" G}] 
    (->> get_map 
         (keep #(when (= (val %) (apply max (vals get_map))) (key %))) 
         rand-nth)))

(defn concensus [file]
  (->> file
       (i/add-derived-column
        :consus
        [:Tun :Aun :Cun :Gun]
        #(get-max  %1 %2 %3 %4))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALL VARIANTS                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def T_matrix (i/dataset [:consus    :Aun  :Tun :Cun  :Gun] 
                         [["A"         0   0.92 1.75  2.00] 
                          ["T"       0.94   0   2.41  1.82] 
                          ["C"       1.81  1.19   0   1.78]
                          ["G"       1.21  1.15 1.81    0 ]]))

(defn poisson [col_var col_val consus c_cov p]
  (let [lambda (->> T_matrix
                    (i/$where {:consus {:$eq consus}})
                    (i/$ col_var)
                    (* (/ c_cov 1000)))]
    (if (< p (st/cdf-poisson col_val :lambda lambda))
      col_val
      0)))
  
(defn pois_correct [col_var col_name p_value file]
  (let [p (- 1 p_value)]
    (->> file
         (i/add-derived-column
          col_name
          [col_var :consus :c_cov]
          #(poisson col_var %1 %2 %3 p)))))


(defn poissonize [p_value file] 
  (->> (pois_correct :Aun :Apois p_value file)
       (pois_correct :Tun :Tpois p_value)
       (pois_correct :Cun :Cpois p_value)
       (pois_correct :Gun :Gpois p_value)))


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

(defn folded-SFS [ref T A C G] 
  (let [f { "T" T "A" A "C" C "G" G}] 
      (second (reverse (sort [ T A C G])))))

(defn SFS [file SFS-type]
  (->> file
       (i/add-derived-column
        :sfs
        [:ref :Tpois :Apois :Cpois :Gpois]
        #(SFS-type  %1 %2 %3 %4 %5))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SUMMARY STATISTICS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn stat-report [file_in]
  (def pied (ii/read-dataset file_in :header true))
  (def nucleotide_diversity (/  (i/sum (i/$ :pie pied)) (i/nrow pied)))
  (def segregation_sites (count (filter #(< 0 %) (i/$ :pie pied))))
  (println "Nucleotide diversity: " nucleotide_diversity)
  (println "Segregating Sites: " segregation_site))
