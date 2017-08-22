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

(defn bm-loc []
  (def L1      (str home "S1.inc"     ))
  (def L10     (str home "S10.inc"    ))
  (def L11     (str home "S11.inc"    ))
  (def L12     (str home "S12.inc"    ))
  (def L13     (str home "S13.inc"    ))
  (def L14     (str home "S14.inc"    ))
  (def L15     (str home "S15.inc"    ))
  (def L16     (str home "S16.inc"    ))
  (def L17     (str home "S17.inc"    ))
  (def L18     (str home "S18.inc"    ))
  (def L23     (str home "S23.inc"    ))
  (def L24     (str home "S24.inc"    ))
  (def L25     (str home "S25.inc"    ))
  (def L26     (str home "S26.inc"    ))
  (def L27     (str home "S27.inc"    ))
  (def L28     (str home "S28.inc"    ))
  (def L29     (str home "S29.inc"    ))
  (def L30     (str home "S30.inc"    ))
  (def L9      (str home "S9.inc"     ))
  (def L05-Pa  (str home "505-Pa.inc" ))
  (def L05-M   (str home "505-M.inc"  ))
  (def L19-Pb  (str home "519-Pb.inc" ))
  (def L19-Pc  (str home "519-Pc.inc" ))
  (def L19-Pd  (str home "519-Pd.inc" ))
  (def L19-S1a (str home "519-S1a.inc"))
  (def L20-Pa  (str home "520-Pa.inc" ))
  (def L20-Pb  (str home "520-Pb.inc" ))
  (def L20-Pc  (str home "520-Pc.inc" ))
  (def L20-S1  (str home "520-S1.inc" ))
  (def L20-S1a (str home "520-S1a.inc"))
  (def L79-Pa  (str home "579-Pa.inc" ))
  (def L79-Pb  (str home "579-Pb.inc" ))
  (def L79-M   (str home "579-M.inc"  ))
  (def L79-S1a (str home "579-S1a.inc"))
  (def L79-S1b (str home "579-S1b.inc")))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-set [file cov]
  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))

(defn bm-sets []
  (bm-loc)
  (def S1      (m-get-set L1  0))
  (def S10     (m-get-set L10 0))

  (def S11     (m-get-set L11 0))
  (def S12     (m-get-set L12 0))
  (def S13     (m-get-set L13 0))
  (def S14     (m-get-set L14 0))

  (def S15     (m-get-set L15 0))
  (def S16     (m-get-set L16 0))
  (def S17     (m-get-set L17 0))
  (def S18     (m-get-set L18 0))
  (def S23     (m-get-set L23 0))

  (def S24     (m-get-set L24 0))
  (def S25     (m-get-set L25 0))
  (def S26     (m-get-set L26 0))
  (def S27     (m-get-set L27 0))
  (def S28     (m-get-set L28 0))
  (def S29     (m-get-set L29 0))
  (def S30     (m-get-set L30 0))
  (def S9      (m-get-set L9  0 ))
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
  (def S79-S1b (m-get-set L79-S1b 0)))


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
                     (map #(add-row (i/$ [:loc :pi :minfr :depth] %)))
                     (map #(i/rename-cols (c-rename %2 %1) %2) (iterate inc 1)))
        enamed  (rest  renamed)]
    (->> renamed
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
              S79-S1a S79-S1b
              S1  S10 S11
              S12 S13
              S14 S15 S16
              S17 S18
              S23 S24 S25
              S26 S27
              S28 S29 S30
              S9])
#_ (def mat
     (bm-loc)
     (bm-sets)
     (pcaM samples))

(defn save-mat [pcaM file-out]
  "/home/yosh/datafiles/incanted_files/SVD15.inc"
  (with-open [f-out (io/writer file-out)]
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




(defn le-filter [file & {:keys [m d]
                          :or   {m 0.01
                                 d 35.0}}]
  (->> file
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  ;;depth4
                   depth5  depth6  depth7
                   depth8  depth9
                  ;; depth10 depth11
                  ;; depth12
                  ;; depth13 depth14
                   depth15 depth16
                   depth17 depth18 depth19
                   depth20 depth21
                   depth22 depth23 depth24
                   depth25 depth26
                   depth27 depth28 depth29
                   depth30 depth31
                   depth32 depth33 depth34
                   depth35]
                  (and (not= nil depth1)
                       (not= nil depth2)  (not= nil depth3)  ;;(not= nil depth4)
                       (not= nil depth5)  (not= nil depth6)  (not= nil depth7)
                       (not= nil depth8)  (not= nil depth9)
                      ;; (not= nil depth10) (not= nil depth11)
                      ;; (not= nil depth12)
                      ;; (not= nil depth13) (not= nil depth14)
                       (not= nil depth15) (not= nil depth16)
                       (not= nil depth17) (not= nil depth18) (not= nil depth19)
                       (not= nil depth20) (not= nil depth21)
                       (not= nil depth22) (not= nil depth23) (not= nil depth24)
                       (not= nil depth25) (not= nil depth26)
                       (not= nil depth27) (not= nil depth28) (not= nil depth29)
                       (not= nil depth30) (not= nil depth31)
                       (not= nil depth32) (not= nil depth33) (not= nil depth34)
                       (not= nil depth35))))
       (i/$where (i/$fn
                  [minfr1
                   minfr2  minfr3  ;;minfr4
                   minfr5  minfr6  minfr7
                   minfr8  minfr9
                   ;; minfr10 minfr11
                   ;; minfr12
                   ;; minfr13 minfr14
                   minfr15 minfr16
                   minfr17 minfr18 minfr19
                   minfr20 minfr21
                   minfr22 minfr23 minfr24
                   minfr25 minfr26
                   minfr27 minfr28 minfr29
                   minfr30 minfr31
                   minfr32 minfr33 minfr34
                   minfr35]
                  (or  (> minfr1 m)
                       (> minfr2 m)  (> minfr3 m) ;; (> minfr4 m)
                       (> minfr5 m)  (> minfr6 m)  (> minfr7 m)
                       (> minfr8 m)  (> minfr9 m)
                       ;;(> minfr10 m) (> minfr11 m)
                       ;;(> minfr12 m)
                       ;;(> minfr13 m) (> minfr14 m)
                       (> minfr15 m) (> minfr16 m)
                       (> minfr17 m) (> minfr18 m) (> minfr19 m)
                       (> minfr20 m) (> minfr21 m)
                       (> minfr22 m) (> minfr23 m) (> minfr24 m)
                       (> minfr25 m) (> minfr26 m)
                       (> minfr27 m) (> minfr28 m) (> minfr29 m)
                       (> minfr30 m) (> minfr31 m)
                       (> minfr32 m) (> minfr33 m) (> minfr34 m)
                       (> minfr35 m))))
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  ;;depth4
                   depth5  depth6  depth7
                   depth8  depth9
                   ;;depth10 depth11
                   ;;depth12
                   ;;depth13 depth14
                   depth15 depth16
                   depth17 depth18 depth19
                   depth20 depth21
                   depth22 depth23 depth24
                   depth25 depth26
                   depth27 depth28 depth29
                   depth30 depth31
                   depth32 depth33 depth34
                   depth35]
                  (and (> depth1  d)
                       (> depth2  d) (> depth3  d) ;;(> depth4  d)
                       (> depth5  d) (> depth6  d) (> depth7  d)
                       (> depth8  d) (> depth9  d)
                       ;;(> depth10 d) (> depth11 d)
                       ;;(> depth12 d)
                       ;;(> depth13 d) (> depth14 d)
                       (> depth15 d) (> depth16 d)
                       (> depth17 d) (> depth18 d) ;;(> depth19 d)
                       (> depth20 d) (> depth21 d)
                       (> depth22 d) (> depth23 d) (> depth24 d)
                       ;;(> depth25 d) (> depth26 d)
                       (> depth27 d) (> depth28 d) (> depth29 d)
                       (> depth30 d) (> depth31 d)
                       ;;(> depth32 d) (> depth33 d) (> depth34 d)
                       ;;(> depth35 d)
                       )))
       (i/$ [:minfr1
             :minfr2  :minfr3  ;;:minfr4
             :minfr5  :minfr6  :minfr7
             :minfr8  :minfr9
             ;;:minfr10 :minfr11
             ;;:minfr12
             ;;:minfr13 :minfr14
             :minfr15 :minfr16
             :minfr17 :minfr18 ;:minfr19
             :minfr20 :minfr21
             :minfr22 :minfr23 :minfr24
             ;;:minfr25 :minfr26
             :minfr27 :minfr28 :minfr29
             :minfr30 :minfr31
             ;;:minfr32 :minfr33 :minfr34
             ;:minfr35
             ])))

 

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
  (let [dataset  (le-filter file)
        old_cols (i/col-names dataset)
        new_cols (vec (map #(keyword (subs (str % "z") 1))
                           old_cols))]
    (->> old_cols
         (reduce #(normalize->mean-zero %2 %1) dataset)
         (i/$ new_cols))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Singular Value Decomposition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn SNP-SVD [file]
  (let [svd  (->> file
                  get-normalized
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

#_(save-mat mat "/mnt/data/datafiles/matrix")

(defn SNP-SVD [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ (range 0 9) 0 projection)
                        (i/$ (range 0 9) 1 projection)
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ (range 9  16) 0 projection)
                      (i/$ (range 9  16) 1 projection))
        (c/add-points (i/$ (range 17 35) 0 projection)
                      (i/$ (range 17 35) 1 projection))
        (i/view))))

