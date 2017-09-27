(ns genome.dreduction
  (require [clojure.java.io   :as io]
           [incanter.core     :as i]
           [incanter.io       :as ii]
           [incanter.stats    :as st]
           [incanter.charts   :as c]
           [clojure.data.csv  :as csv]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;PCA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn  c-rename [file ser_num]
  "This one is for adding a serial number at the end of a col name"
  (let [dont-change #{ :loc}
        old-cols (apply vector (remove #(contains? dont-change %)
                                       (i/col-names file)))
        new-cols (->> old-cols
                      (map #(keyword (subs (str % ser_num) 1)))
                      (apply vector))]
    (->> new-cols
         (interleave old-cols)
         (apply assoc {}))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Reduce

(defn PCA-matrix [samples]
  (let [renamed (->> samples
                     (map #(add-row (i/$ [:loc :pi :mf :depth] %)))
                     (map #(i/rename-cols (c-rename %2 %1) %2) (iterate inc 1)))
        enamed  (rest  renamed)]
    (->> renamed
         (reduce #(i/$join [:loc :loc] %1 %2) (map enamed))
         (i/$where {:pi1 {:$ne nil}}))))

(defn zero-fr [m file]
  "create a new collumnt for one dataset with freq > m"
  (->> file
       (i/add-derived-column
        :mf
        [:minfr]
        #(if (> m  %) 0.0 %))))

(defn get-zerowed-fr [m samples]
  (PCA-matrix (map #(zero-fr m %) samples)))

#_(def mat (get-zerowed-fr 0.05 samples))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Dimention reduction analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn clean-low [file]
  "cleans low"
  (->> file
       (i/$where (i/$fn [pi] (> pi 0.0)))
       (i/add-derived-column
        :var-num
        [:cov :minfr]
        #(* %1 %2))))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Normalizing covariance matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn  normalize->mean-zero [col_name file]
  "snp sould be the string of the keyword"
  (let [m (st/mean (i/$ col_name file))]
    (->> file
         (i/add-derived-column
          (keyword (subs (str col_name "z") 1))
          [col_name]
          #(- %1 m)))))

(defn get-normalized [file]
  (let [old_cols (i/col-names file)
        new_cols (vec (map #(keyword (subs (str % "z") 1))
                           old_cols))]
    (->> old_cols
         (reduce #(normalize->mean-zero %2 %1) file)
         (i/$ new_cols))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Singular Value Decomposition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn SNP-SVD [file]
  (let [svd  (->> file
                  i/to-matrix
                  i/trans
                  i/decomp-svd)
        dims 2
        u    (i/$
             (range dims) (:U svd))
        s    (i/diag (take dims
                           (:S svd)))
        v   (i/trans (i/$ (range dims) (:V svd)))]
    (i/mmult u s)))

(defn save-mat [file file-out]
  "/home/yosh/datafiles/incanted_files/SVDs.inc"
  (let [SVDs (SNP-SVD file)]
    (with-open [f-out (io/writer file-out)]
      (csv/write-csv f-out (i/to-list SVDs)))))



