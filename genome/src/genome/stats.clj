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

;Creates a sliding window from a column :pi in file and adds a new col :pislide 
(defn pi-slide [file column win_size]
  (let [winset (ii/read-dataset file :header true)]
    (i/add-column
     :sliding
     (->> (i/$ column winset)
          (partition win_size 1)
          (map #(/ (apply + %) win_size))
          (concat (take (dec win_size) (repeat 0))))
     winset)))

(defn pi-slide-try [file column win_size]
  (let [winset (ii/read-dataset file :header true)]
    (i/add-column
     :slidingTry
     (->> (i/$ column winset)
          (z/roll-mean win_size)
          (concat (take (dec win_size) (repeat 0))))
     winset)))


;in the future I would like to move all stats here from database
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

(defn poisson [col_var c_cov err_mil p]
  (let [lambda (* (/ c_cov 1000) err_mil)] ;error corrected to lambda by coverage
  (if (> p (st/cdf-poisson col_var :lambda lambda))
    col_var
    0))

  
(defn pois_correct [file col_var col_name err_mil p_value]
  (let [p (- 1 p_value)]
    (->> file
         (i/add-derived-column
          col_name
          [col_var :c_cov]
          #(poisson %1 %2 err_mil p)))))

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
  (let [winset (ii/read-dataset file :header true)]
    (->> winset
         (i/add-derived-column
          :sfs
          [:ref :Tun :Aun :Cun :Gun]
          #(SFS-type  %1 %2 %3 %4 %5)))))

(defn stat-report [file_in]
  (def pied (ii/read-dataset file_in :header true))
  (def nucleotide_diversity (/  (i/sum (i/$ :pie pied)) (i/nrow pied)))
  (def segregation_sites (count (filter #(< 0 %) (i/$ :pie pied))))
  (println "Nucleotide diversity: " nucleotide_diversity)
  (println "Segregating Sites: " segregation_site))
