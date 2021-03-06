(ns genome.pop
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

(defn get-pi [cov T A C G] 
"Calculates pi for each row"
  (if (>=  cov 2) 
    (double (/(+ (* T A) (* T C)
                 (* T G) (* A C) 
                 (* A G) (* G C))
              (/ (* cov (- cov 1))
                 2)))
    0.0))

(def pi_un [:pi_un [:cov_un :Tun :Aun :Cun :Gun]]);for variants calling
(def pi    [:pi    [:depth  :T   :A   :C   :G  ]]);after variants

(defn pise [pi_type file]
  (->> file
       (i/add-derived-column
         (first pi_type)
         (last pi_type) 
         #(get-pi %1 %2 %3 %4 %5))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SLIDING WINDOW
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn slide-mean [file scanned_column new_column win_size]
"Creates a sliding window from a column :pi in file and adds a new col :pislide"
  (i/add-column
   new_column
   (->> (i/$ scanned_column file)
        (partition win_size 1)
        (map #(/ (apply + %) win_size))
        (concat (take (dec win_size) (repeat 0))))
   file))
(def m-slide-mean (memoize slide-mean))


(defn glide-mean [file scanned_column new_column win_size]
"Very similar to 'slide-mean' with a built in function"
  (i/add-column
   new_column
   (->> (i/$ scanned_column file)
        (z/roll-apply win_size) ;instead of roll-min we can use roll-apply func
        (concat (take (dec win_size) (repeat 0))))
   file))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ALLELE FREQUENCY SPECTRA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn unfolded-SFS [ref T A C G] 
"Calculates unfolded site allele frequency for prevalent minor alleles"
  (if-not (= ref "-")
    (let [f { "T" T "A" A "C" C "G" G}] 
      (apply max (filter #(not= (f ref) %) [ T A C G])))
    "-")) 


(defn multi-SFS [ref T A C G] 
"Calculates unfolded site allele frequency for all minor alleles"
  (if-not (= ref "-")
    (let [f { "T" T "A" A "C" C "G" G}] 
      (i/sum (filter #(not= (f ref) %) [ T A C G])))
    "-"))

(defn folded-SFS [mean_cov ref cov_p T A C G] 
"Calculates folded site allele frequency"
  (if (= 0.0 cov_p)
    0.0
    (->> [T A C G]
         sort
         reverse
         second
         (* (/ mean_cov cov_p))
         double)))

(defn SFS [SFS-type file]
  (let [mean_cov (st/mean (i/$ :depth file))]
    (->> file
         (i/add-derived-column
          :sfs
          [:maj_un+ :depth :T :A :C :G]
          #(SFS-type mean_cov %1 %2 %3 %4 %5 %6)))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;BINNING
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn bin-sfs [n-bins file] ;bin range 0-14 into 5 bins (bin 5 (range 15))
  (let [sfs_0      (i/nrow   (i/$where (i/$fn [sfs] (= 0.0 sfs)) file))
        bin-array  (i/$ :sfs (i/$where (i/$fn [sfs] (< 0.0 sfs)) file))
        min-freq   (apply min bin-array)
        max-freq   (apply max bin-array)
        range-freq (- max-freq min-freq)
        bin-fn     (fn [freq](-> freq
                                 (- min-freq)
                                 (/ range-freq)
                                 (* n-bins)
                                 (int)
                                 (inc)
                                 (min n-bins)))]
    (conj (->> (map bin-fn bin-array)
               frequencies
               sort)
          [0 sfs_0])))

