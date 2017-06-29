(ns genome.dreduction
  (require [clojure.java.io   :as io]
           [incanter.core     :as i]
           [incanter.io       :as ii]
           [incanter.stats    :as st]
           [incanter.charts   :as c]
           [clojure.data.csv  :as csv]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;
;; For PC
;;(def home "/home/yosh/datafiles/incanted_files/")

;;;;;;;;;;;;;;;;;;;;;;;
;; For Server
(def home "/mnt/data/datafiles/incanted_files/")


(def L05-Pa  (str home "505-Pa.inc"))
(def L05-M   (str home "505-M.inc"))

(def L19-Pb  (str home "519-Pb.inc"))
(def L19-Pc  (str home "519-Pc.inc"))
(def L19-Pd  (str home "519-Pd.inc"))
(def L19-S1a (str home "519-S1a.inc"))

(def L20-Pa  (str home "520-Pa.inc"))
(def L20-Pb  (str home "520-Pb.inc"))
(def L20-Pc  (str home "520-Pc.inc"))
(def L20-S1  (str home "520-S1.inc"))
(def L20-S1a (str home "520-S1a.inc"))

(def L79-Pa  (str home "579-Pa.inc"))
(def L79-Pb  (str home "579-Pb.inc"))
(def L79-M   (str home "579-M.inc"))
(def L79-S1a (str home "579-S1a.inc"))
(def L79-S1b (str home "579-S1b.inc"))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-set [file cov]
  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))

(defn les-sets [])
(def S05-Pa  (m-get-set L05-Pa  0))
(def S05-M   (m-get-set L05-M   0))

(def S19-Pb  (m-get-set L19-Pb  0))
(def S19-Pc  (m-get-set L19-Pc  0))
(def S19-Pd  (m-get-set L19-Pd  0))
(def S19-S1a (m-get-set L19-S1a 0))

(def S20-Pa  (m-get-set L20-Pa  0))
(def S20-Pb  (m-get-set L20-Pb  0))
(def S20-Pc  (m-get-set L20-Pc  0))
(def S20-S1a (m-get-set L20-S1a 0))
(def S20-S1b (m-get-set L20-S1  0))

(def S79-Pa  (m-get-set L79-Pa  0))
(def S79-Pb  (m-get-set L79-Pb  0))
(def S79-M   (m-get-set L79-M   0))
(def S79-S1a (m-get-set L79-S1a 0))
(def S79-S1b (m-get-set L79-S1b 0))


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
                     (map #(add-row (i/$ [:loc :pi :minfr :depth] %)))
                     (map #(i/rename-cols (c-rename %2 %1) %2) (iterate inc 1)))
        renamed (first renamed)
        enamed  (rest  renamed)]
    (->> r
         (reduce #(i/$join [:loc :loc] %1 %2) (map enamed))
         (i/$where {:pi1 {:$ne nil}}))))


(def pcaM    (memoize PCA-matrix))
(def samples [S05-Pa
              S19-Pb  S19-Pc  S19-Pd
              S20-Pa  S20-Pb  S20-Pc
              S79-Pa  S79-Pb
              S05-M   S79-M
              S19-S1a
              S20-S1a S20-S1b
              S79-S1a S79-S1b])


(defn save-mat [pcaM file-out]
  "/home/yosh/datafiles/incanted_files/SVD15.inc"
  (with-open [f-out (io/writer file_out)]
    (csv/write-csv f-out [(map name (i/col-names pcaM))])
    (csv/write-csv f-out (i/to-list pcaM))))


(defn get-mat [file]
  "open an csv.inc file"
  (ii/read-dataset file :header true))
(def m-get-set (memoize get-set))


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



(defn le-filter [file]
  (->> file
       (i/$where (i/$fn [depth1  depth2  depth3  depth4
                         depth5  depth6  depth7  depth8
                         depth9  depth10 depth11 depth12
                         depth13 depth14 depth15 depth16]
                        (and (not= nil depth1)  (not= nil depth2)
                             (not= nil depth3)  (not= nil depth4)
                             (not= nil depth5)  (not= nil depth6)
                             (not= nil depth7)  (not= nil depth8)
                             (not= nil depth9)  (not= nil depth10)
                             (not= nil depth11) (not= nil depth12)
                             (not= nil depth13) (not= nil depth14)
                             (not= nil depth15) (not= nil depth16))))
       (i/$where (i/$fn [minfr1  minfr2  minfr3  minfr4
                         minfr5  minfr6  minfr7  minfr8
                         minfr9  minfr10 minfr11 minfr12
                         minfr13 minfr14 minfr15 minfr16]
                        (or  (> minfr1 0.0035)  (> minfr2 0.0035)
                             (> minfr3 0.0035)  (> minfr4 0.0035)
                             (> minfr5 0.0035)  (> minfr6 0.0035)
                             (> minfr7 0.0035)  (> minfr8 0.0035)
                             (> minfr9 0.0035)  (> minfr10 0.0035)
                             (> minfr11 0.0035) (> minfr12 0.0035)
                             (> minfr13 0.0035) (> minfr14 0.0035)
                             (> minfr15 0.0035) (> minfr16 0.0035))))
       (i/$where (i/$fn [depth1  depth2  depth3  depth4
                         depth5  depth6  depth7  depth8
                         depth9  depth10 depth11 depth12
                         depth13 depth14 depth15 depth16]
                        (and (> depth1  20.0) (> depth2  20.0)
                             (> depth3  20.0) (> depth4  20.0)
                             (> depth5  20.0) (> depth6  20.0)
                             (> depth7  20.0) (> depth8  20.0)
                             (> depth9  20.0) (> depth10 20.0)
                             (> depth11 20.0) (> depth12 20.0)
                             (> depth13 20.0) (> depth14 20.0)
                             (> depth15 20.0) (> depth16 20.0))))
       (i/$ [:minfr1  :minfr2  :minfr3  :minfr4
             :minfr5  :minfr6  :minfr7  :minfr8
             :minfr9  :minfr10 :minfr11 :minfr12
             :minfr13 :minfr14 :minfr15 :minfr16])))

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

(defn get-normalized [from1 to1 file]
  (let [dataset  (i/$ (range from1 to1) :all (le-filter file))
        old_cols (i/col-names dataset)
        new_cols (vec (map #(keyword (subs (str % "z") 1))
                           old_cols))]
    (->> old_cols
         (reduce #(normalize->mean-zero %2 %1) dataset)
         (i/$ new_cols))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Singular Value Decomposition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn SNP-SVD [file from1 to1]
  (let [svd (->> file
                 (get-normalized from1 to1)
                 i/to-matrix
                 i/trans
                 i/decomp-svd)
        dims 2
        u (i/$
           (range dims) (:U svd))
        s (i/diag (take dims
                        (:S svd)))
        v (i/trans (i/$ (range dims) (:V svd)))
        projection (i/mmult u s)]
    (-> (c/scatter-plot (i/$ (range 0 9) 0 projection)
                        (i/$ (range 0 9) 1 projection)
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ (range 9 16) 0 projection)
                      (i/$ (range 9 16) 1 projection))
        (i/view))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Pricipal Components Analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn SNPxs [file]
  (let [data       (->> file
                        le-filter
                        i/to-matrix
                        i/trans)
        components (st/principal-components data)
        pc1        (i/$ 0 (:rotation components))]
    (i/mmult data pc1)))

(defn SNPys [file]
  (let [data       (->> file
                        le-filter
                        i/to-matrix
                        i/trans)
        components (st/principal-components data)
        pc2        (i/$ 1 (:rotation components))]
    (i/mmult data pc2)))


(defn SNPCA [xs ys]
  (-> (c/scatter-plot (i/$ (range 0 9) 0 xs)
                      (i/$ (range 0 9) 0 ys)
                      :x-label "Principle Component 1"
                      :y-label "Principle Component 2")
      (c/add-points (i/$  (range 9 10) 0 xs)
                    (i/$  (range 9 10) 0 ys))
      (c/add-points (i/$  (range 10 11) 0 xs)
                    (i/$  (range 10 11) 0 ys))
      (c/add-points (i/$  (range 11 13) 0 xs)
                    (i/$  (range 11 13) 0 ys))
      (c/add-points (i/$  (range 13 14) 0 xs)
                    (i/$  (range 13 14) 0 ys))
      (i/view)))

(def h (PCA-matrix [(i/$ (range 0 1000) :all S19-Pb)
                    (i/$ (range 0 1000) :all S19-Pc)
                    (i/$ (range 0 1000) :all S19-Pc)]))


