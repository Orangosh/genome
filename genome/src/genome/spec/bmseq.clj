(ns genome.spec.all
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [incanter.charts   :as c  ]
           [genome.dna2aa     :as da ]
           [genome.stats      :as gs ]
           [genome.consvar    :as cv ]
           [genome.pop        :as p  ]
           [clojure.data.csv  :as csv]
           [genome.view       :as v  ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;
;; For PC
;;(def home "/home/yosh/data")

;;;;;;;;;;;;;;;;;;;;;;;
;; For Server
(def home "/mnt/data/")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn all-cmv-loc []
  (def L05-Pa  (str home "hcmv/incanted_files/505-Pa.inc" ))
  (def L05-M   (str home "hcmv/incanted_files/505-M.inc"  ))

  (def L19-Pb  (str home "hcmv/incanted_files/519-Pb.inc" ))
  (def L19-Pc  (str home "hcmv/incanted_files/519-Pc.inc" ))
  (def L19-Pd  (str home "hcmv/incanted_files/519-Pd.inc" ))
  (def L19-S1a (str home "hcmv/incanted_files/519-S1a.inc"))

  (def L20-Pa  (str home "hcmv/incanted_files/520-Pa.inc" ))
  (def L20-Pb  (str home "hcmv/incanted_files/520-Pb.inc" ))
  (def L20-Pc  (str home "hcmv/incanted_files/520-Pc.inc" ))
  (def L20-S1  (str home "hcmv/incanted_files/520-S1.inc" ))
  (def L20-S1a (str home "hcmv/incanted_files/520-S1a.inc"))

  (def L79-Pa  (str home "hcmv/incanted_files/579-Pa.inc" ))
  (def L79-Pb  (str home "hcmv/incanted_files/579-Pb.inc" ))
  (def L79-M   (str home "hcmv/incanted_files/579-M.inc"  ))
  (def L79-S1a (str home "hcmv/incanted_files/579-S1a.inc"))
  (def L79-S1b (str home "hcmv/incanted_files/579-S1b.inc"))

  ;;bm-sets
  (def L1      (str home "bmseq/incanted_files/S1.inc"     ))

  (def L10     (str home "bmseq/incanted_files/S10.inc"    ))
  (def L11     (str home "bmseq/incanted_files/S11.inc"    ))
  (def L12     (str home "bmseq/incanted_files/S12.inc"    ))
  (def L13     (str home "bmseq/incanted_files/S13.inc"    ))
  (def L14     (str home "bmseq/incanted_files/S14.inc"    ))
  (def L15     (str home "bmseq/incanted_files/S15.inc"    ))
  (def L16     (str home "bmseq/incanted_files/S16.inc"    ))
  (def L17     (str home "bmseq/incanted_files/S17.inc"    ))
  (def L18     (str home "bmseq/incanted_files/S18.inc"    ))

  (def L23     (str home "bmseq/incanted_files/S23.inc"    ))
  (def L24     (str home "bmseq/incanted_files/S24.inc"    ))
  (def L25     (str home "bmseq/incanted_files/S25.inc"    ))
  (def L26     (str home "bmseq/incanted_files/S26.inc"    ))
  (def L27     (str home "bmseq/incanted_files/S27.inc"    ))
  (def L28     (str home "bmseq/incanted_files/S28.inc"    ))
  (def L29     (str home "bmseq/incanted_files/S29.inc"    ))
  (def L30     (str home "bmseq/incanted_files/S30.inc"    ))

  (def L9      (str home "bmseq/incanted_files/S9.inc"     )))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-set [file cov]
  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))


(defn all-cmv-sets []
  (bmseq-loc)
  (def samples  {:S1      (m-get-set L1  0)
                 :S10     (m-get-set L10 0)

                 :S11     (m-get-set L11 0)
                 :S12     (m-get-set L12 0)
                 :S13     (m-get-set L13 0)
                 :S14     (m-get-set L14 0)

                 :S15     (m-get-set L15 0)
                 :S16     (m-get-set L16 0)
                 :S17     (m-get-set L17 0)
                 :S18     (m-get-set L18 0)
                 :S23     (m-get-set L23 0)

                 :S24     (m-get-set L24 0)
                 :S25     (m-get-set L25 0)
                 :S26     (m-get-set L26 0)
                 :S27     (m-get-set L27 0)
                 :S28     (m-get-set L28 0)
                 :S29     (m-get-set L29 0)
                 :S30     (m-get-set L30 0)
                 :S9      (m-get-set L9  0 )
                 :S05-Pa  (m-get-set L05-Pa  0)
                 :S05-M   (m-get-set L05-M   0)

                 :S19-Pb  (m-get-set L19-Pb  0)
                 :S19-Pc  (m-get-set L19-Pc  0)
                 :S19-Pd  (m-get-set L19-Pd  0)
                 :S19-S1a (m-get-set L19-S1a 0)

                 :S20-Pa  (m-get-set L20-Pa  0)
                 :S20-Pb  (m-get-set L20-Pb  0)
                 :S20-Pc  (m-get-set L20-Pc  0)
                 :S20-S1a (m-get-set L20-S1a 0)
                 :S20-S1b (m-get-set L20-S1  0)

                 :S79-Pa  (m-get-set L79-Pa  0)
                 :S79-Pb  (m-get-set L79-Pb  0)
                 :S79-M   (m-get-set L79-M   0)
                 :S79-S1a (m-get-set L79-S1a 0)
                 :S79-S1b (m-get-set L79-S1b 0)}))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;FOR PCA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def mat
  (genome.dreduction/pcaM samples))


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

#_(genome.dreduction/save-mat (le-filter mat) "/mnt/data/bmseq/all-cmv-matrix")

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



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET WINDOWED SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#_(defn win-100 [file]
    (p/m-slide-mean file :pi :pi_slide 100))

(defn win-100 [file]
  (p/m-slide-mean file :minfr :minfr_slide 100))

(defn all-cmv-win []
  (def W1      (win-100 (m-get-set L1      20)))
  (def W10     (win-100 (m-get-set L10     20)))

  (def W11     (win-100 (m-get-set L11     20)))
  (def W12     (win-100 (m-get-set L12     20)))
  (def W13     (win-100 (m-get-set L13     20)))
  (def W14     (win-100 (m-get-set L14     20)))

  (def W15     (win-100 (m-get-set L15     20)))
  (def W16     (win-100 (m-get-set L16     20)))
  (def W17     (win-100 (m-get-set L17     20)))
  (def W18     (win-100 (m-get-set L18     20)))
  (def W23     (win-100 (m-get-set L23     20)))

  (def W24     (win-100 (m-get-set L24     20)))
  (def W25     (win-100 (m-get-set L25     20)))
  (def W26     (win-100 (m-get-set L26     20)))
  (def W27     (win-100 (m-get-set L27     20)))
  (def W28     (win-100 (m-get-set L28     20)))
  (def W29     (win-100 (m-get-set L29     20)))
  (def W30     (win-100 (m-get-set L30     20)))
  (def W9      (win-100 (m-get-set L9      20)))
  (def W05-Pa  (win-100 (m-get-set L05-Pa  20)))
  (def W05-M   (win-100 (m-get-set L05-M   20)))

  (def W19-Pb  (win-100 (m-get-set L19-Pb  20)))
  (def W19-Pc  (win-100 (m-get-set L19-Pc  20)))
  (def W19-Pd  (win-100 (m-get-set L19-Pd  20)))
  (def W19-S1a (win-100 (m-get-set L19-S1a 20)))

  (def W20-Pa  (win-100 (m-get-set L20-Pa  20)))
  (def W20-Pb  (win-100 (m-get-set L20-Pb  20)))
  (def W20-Pc  (win-100 (m-get-set L20-Pc  20)))
  (def W20-S1a (win-100 (m-get-set L20-S1a 20)))
  (def W20-S1  (win-100 (m-get-set L20-S1  20)))

  (def W79-Pa  (win-100 (m-get-set L79-Pa  20)))
  (def W79-Pb  (win-100 (m-get-set L79-Pb  20)))
  (def W79-M   (win-100 (m-get-set L79-M   20)))
  (def W79-S1a (win-100 (m-get-set L79-S1a 20)))
  (def W79-S1b (win-100 (m-get-set L79-S1b 20))))

(defn look
  "create a graph of the value comparison"
  ([column file]
   (i/view (c/xy-plot
            :loc
            column
            :x-label "frequency"
            :y-label "Location"
            :title   "Change in nucelotide diversity per site between 2 samples"
;:legend true
            :data file)))

  ([column early>inter early>late]
   (-> (c/xy-plot
        :loc
        column
        :x-label "frequency"
        :y-label "Location"
        :title   "Change in nucelotide diversity per site between 2 samples"
        :data early>inter)
       (c/add-lines
        :loc
        column
        :data early>late)
       (i/view)))

  ([column early>inter early>late early>mom_sib]
    (-> (c/xy-plot
         :loc
         column
         :x-label "frequency"
         :y-label "Location"
         :title   "Change in nucelotide diversity per site between 2 samples"
         :data early>inter)
        (c/add-lines
         :loc
         column
         :data early>late)
        (c/add-lines
         :loc
         column
         :data early>mom_sib)
        (i/view)))

  ([column file1 file2 file3 file4 file5 file6]
   (-> (c/xy-plot
        :loc
        column
        :x-label "frequency"
        :y-label "Location"
        :title   "Change in nucelotide diversity per site between 2 samples"
        :data file1)
       (c/add-lines
        :loc
        column
        :data file2)
       (c/add-lines
        :loc
        column
        :data file3)
       (c/add-lines
        :loc
        column
        :data file4)
       (c/add-lines
        :loc
        column
        :data file5)
       (c/add-lines
        :loc
        column
        :data file6)
       (i/view))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Gneral description
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Gneral descriptive functions


(defn cutoff [file]
  (i/nrow (i/$where (i/$fn [minfr pi]
                           (and (> minfr 0.0035)
                                (> pi 0.0))) file)))
(defn sum-cov [file]
  (i/sum (i/$ :depth file)))

(defn sum-pi [file]
  (i/sum (i/$ :pi (i/$where (i/$fn [minfr] (> minfr 0.0035)) file))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Samples dataset
(defn samples []
  (->> (i/dataset
        [:sample  :player   :time-pt]
        [[S1  "Primary" 1]
         [S10 "Mother"  1]
         [S11 "Primary" 1]
         [S12 "Primary" 2]
         [S13 "Primary" 3]
         [S14 "Sibling" 1]
         [S15 "Primary" 1]
         [S16 "Primary" 2]
         [S17 "Primary" 3]
         [S18 "Sibling" 1]
         [S23 "Sibling" 2]
         [S24 "Primary" 1]
         [S25 "Primary" 2]
         [S26 "Mother"  0]
         [S27 "Sibling" 1]
         [S28 "Sibling" 2]
         [S29 "Sibling" 2]
         [S30 "Sibling" 2]
         [S9  "Sibling" 2] ])
       (i/add-derived-column
        :name
        [:sample]
        #(subs (first (i/$ :r_seq  %)) 3))
       (i/add-derived-column
        :mean-cov
        [:sample]
        #(/ (sum-cov %)
            (i/nrow  %)))
       (i/add-derived-column
        :cov>20
        [:sample]
        #(i/$where (i/$fn [cov] (> cov 20)) %))
       (i/add-derived-column
        :n-seg
        [:cov>20]
        #(if (> (i/nrow %) 0)
           (/ (double (cutoff %)) (i/nrow %)) 0.0))
       (i/add-derived-column
        :nuc-div
        [:cov>20]
        #(if (> (i/nrow %) 0)
           (/ (sum-pi %) (i/nrow %)) 0.0))
       (i/add-derived-column
        :sfs
        [:cov>20]
        #(if (> (i/nrow %) 0) (p/bin-sfs 10 (da/get-synonymous %)) 0))
       (i/$ [:name :player :time-pt :mean-cov :n-seg :nuc-div :sample :cov>20 :sfs])))
(def p-samples (memoize samples))
