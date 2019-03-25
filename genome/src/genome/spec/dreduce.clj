(ns genome.spec.dreduce
  (require [clojure.java.io   :as io]
           [incanter.core     :as i ]
           [incanter.io       :as ii]
           [incanter.stats    :as st]
           [incanter.charts   :as c ]
           [clojure.data.csv  :as csv]
           [genome.spec.getseqs :refer :all]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Replace Nucleotides with numbers create a matrix and centralizes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn str>num
  [x]
  (case x
    \- 0
    \A 1 \a 1
    \C 2 \c 2
    \T 3 \t 3
    \G 4 \g 4
    0))

(defn str-seq>num-seq
  "changes map value form a nucleotide stringe to a seq of numbers"
  [[k v]]
  (mapv #(str>num %) v))

(defn get-matrix
  "creats a matrix without refernce column and reneams columes according to samples"
  [selected sample_name]
  (let [merlinless (dissoc selected ">Merlin")]
    (->> (apply i/conj-cols
               (map #(str-seq>num-seq %) merlinless))
          (#(i/rename-cols (apply hash-map (interleave (i/col-names %) (keys merlinless)))%)))))



(defn add-mean-col [file]
    "calculates mean for each dimantion"
  (->> file
       (i/add-derived-column
        :mean
        (i/col-names file)
        #(double (/ (+ %1  %2  %3  %4
                       %5  %6  %7  %8
                       %9  %10 %11 %12
                       %13 %14 %15 %16)  16)))))

(defn normalize-a-col
  "reduces mean form a sample"
  [file col]
  (->> file
       (i/add-derived-column
        (keyword col)
        [col :mean]
        #(- %1 %2 ))))                     

(defn normalize
  "calculates mean for each dimantiona and reduce it from all members"
  [file]
  (->> (reduce #(normalize-a-col %1 %2) file (i/col-names file))
       (i/$ [:>505-Pa
             :>519-Pb :>519-Pc :>519-Pd
             :>520-Pa :>520-Pb :>520-Pc
             :>579-Pb :>579-Pa
             :>505-M 
             :>519-S1a 
             :>520-S1a :>520-S1b
             :>579-M :>579-S1a :>579-S1b ])))

(def cmv-normalized
  (->> (get-matrix hcmv ">Merlin")
       (add-mean-col)
       (normalize)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Dimention reduction analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;PCA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn PCA [file]
  (let [mat  (->> file
                 i/to-matrix)
        covm (i/mult mat (i/trans mat))]
    (st/principal-components covm)))

(def SVD (snp-SVD cmv-normalized))

(defn save-mat [file file-out]
  "/mnt/data/hcmv/consensus_analysis/pca/cmv-pca"
  (let [cmv-pca    (SNP-SVD file)
        components (cmv-pca :rotation)
        pc1 (i/sel components :cols 0)
        pc2 (i/sel components :cols 1)
        x2  (i/mmult (i/to-matrix cmv-normalized) pc2)
        x2  (i/mmult (i/to-matrix cmv-normalized) pc2)]
    (with-open [f-out (io/writer file-out)]
      (csv/write-csv f-out (i/to-list cmv-pca)))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Singular Value Decomposition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn SNP-SVD [file]
  (let [svd  (->> file
                  i/to-matrix
                  #_i/trans
                  i/decomp-svd)
        dims 2
        u    (i/$
             (range dims) (:U svd))
        s    (i/diag (take dims
                           (:S svd)))
        v   (i/trans (i/$ (range dims) (:V svd)))]
    (i/mmult u s)))


