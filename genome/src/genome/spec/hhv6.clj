(ns genome.spec.hhv6
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [incanter.charts   :as c  ]
           [genome.dna2aa     :as da ]
           [genome.stats      :as gs ]

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


(defn hhv6-loc []
  ;;hhv6-sets 
  (def L37-Pa  (str home "hhv6/incanted_files/537-Pa.inc" ))
  (def L37-S2b (str home "hhv6/incanted_files/537-S2b.inc"))        
  (def L37-Pb  (str home "hhv6/incanted_files/537-Pb.inc" ))
  (def L37-S3  (str home "hhv6/incanted_files/537-S3.inc" ))       
  (def L37-Pc  (str home "hhv6/incanted_files/537-Pc.inc" ))            
  (def L37-S2a (str home "hhv6/incanted_files/537-S2a.inc"))

  (def L42-Mb  (str home "hhv6/incanted_files/542-Mb.inc" ))
  (def L42-Ma  (str home "hhv6/incanted_files/542-Ma.inc" ))
  (def L42-Pa  (str home "hhv6/incanted_files/542-Pa.inc" ))
  (def L42-Pb  (str home "hhv6/incanted_files/542-Pb.inc" ))
  (def L42-Pc  (str home "hhv6/incanted_files/542-Pc.inc" ))
  (def L42-S1a (str home "hhv6/incanted_files/542-S1a.inc"))
  (def L42-S1b (str home "hhv6/incanted_files/542-S1b.inc"))
  (def L42-S1c (str home "hhv6/incanted_files/542-S1c.inc"))

  (def L43-Pb  (str home "hhv6/incanted_files/543-Pb.inc" ))
  (def L43-S1a (str home "hhv6/incanted_files/543-S1a.inc"))  
  (def L43-S1b (str home "hhv6/incanted_files/543-S1b.inc"))

  (def L72-M   (str home "hhv6/incanted_files/572-M.inc"  ))
  (def L72-Pa  (str home "hhv6/incanted_files/572-Pa.inc" )))


(defn get-set [file cov]

  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))


(defn hhv6-sets []
  (hhv6-loc)
  (def samples {:H37-Pa  (m-get-set L37-Pa  0)
                :H37-S2b (m-get-set L37-S2b 0)        
                :H37-Pb  (m-get-set L37-Pb  0)
                :H37-S3  (m-get-set L37-S3  0)       
                :H37-Pc  (m-get-set L37-Pc  0)            
                :H37-S2a (m-get-set L37-S2a 0)

                :H42-Mb  (m-get-set L42-Mb  0)
                :H42-Ma  (m-get-set L42-Ma  0)
                :H42-Pa  (m-get-set L42-Pa  0)
                :H42-Pb  (m-get-set L42-Pb  0)
                :H42-Pc  (m-get-set L42-Pc  0)
                :H42-S1a (m-get-set L42-S1a 0)
                :H42-S1b (m-get-set L42-S1b 0)
                :H42-S1c (m-get-set L42-S1c 0)

                :H43-Pb  (m-get-set L43-Pb  0)
                :H43-S1a (m-get-set L43-S1a 0)  
                :H43-S1b (m-get-set L43-S1b 0)

                :H72-M   (m-get-set L72-M   0)
                :H72-Pa  (m-get-set L72-Pa  0)}))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;


(defn le-filter [file & {:keys [m d]
                          :or   {m 0.01
                                 d 13.0}}]
  (->> file
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  depth4
                   depth5  depth6  depth7
                   depth8  depth9
                   depth10 depth11
                   depth12 depth13
                   depth14  
                   depth15 depth16
                   depth17 depth18 depth19] 
                  (and (not= nil depth1)
                       (not= nil depth2)  (not= nil depth3)  (not= nil depth4)
                       (not= nil depth5)  (not= nil depth6)  (not= nil depth7)
                       (not= nil depth8)  (not= nil depth9)
                       (not= nil depth10) (not= nil depth11)
                       (not= nil depth12)
                       (not= nil depth13) (not= nil depth14)
                       (not= nil depth15) (not= nil depth16)
                       (not= nil depth17) (not= nil depth18) (not= nil depth19)))) 
       (i/$where (i/$fn
                  [minfr1
                   minfr2  minfr3  minfr4

                   minfr5  minfr6  minfr7
                   minfr8  minfr9
                   minfr10 minfr11
                   minfr12
                   minfr13 minfr14
                   minfr15 minfr16
                   minfr17 minfr18 minfr19]
                  (or  (> minfr1 m)
                       (> minfr2 m)  (> minfr3 m)  (> minfr4 m)
                       (> minfr5 m)  (> minfr6 m)  (> minfr7 m)
                       (> minfr8 m)  (> minfr9 m)
                       (> minfr10 m) (> minfr11 m)
                       (> minfr12 m)
                       (> minfr13 m) (> minfr14 m)
                       (> minfr15 m) (> minfr16 m)
                       (> minfr17 m) (> minfr18 m) (> minfr19 m)))) 
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  depth4
                   depth5  depth6  depth7
                   depth8  depth9
                   depth10 depth11
                   depth12
                   depth13 depth14
                   depth15 depth16
                   depth17 depth18 depth19]
                  (and (> depth1  d)
                       (> depth2  d) (> depth3  d) (> depth4  d)
                       (> depth5  d) (> depth6  d) (> depth7  d)
                       (> depth8  d) (> depth9  d)
                       (> depth10 d) (> depth11 d)
                       (> depth12 d)
                       (> depth13 d) (> depth14 d)
                       (> depth15 d) (> depth16 d)
                       (> depth17 d) (> depth18 d) (> depth19 d))))
       (i/$ [:minfr1
             :minfr2  :minfr3  :minfr4
             :minfr5  :minfr6  :minfr7
             :minfr8  :minfr9
             :minfr10 :minfr11
             :minfr12
             :minfr13 :minfr14
             :minfr15 :minfr16
             :minfr17 :minfr18 :minfr19])))

(defn HHV6-SVD-primary [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 1 2 8 9 10 14 18] 0 projection)
                        (i/$ [0 1 2 8 9 10 14 18] 1 projection)
                        :title "HHV6 primary vs mother vs sibling"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [3 4 5 11 12 13 15 16]  0 projection)
                      (i/$ [3 4 5 11 12 13 15 16]  1 projection))
        (c/add-points (i/$ [6 7 17] 0 projection)
                      (i/$ [6 7 17] 1 projection))
        (i/view))))


(defn HHV6-SVD-family [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 1 2 3 4 5] 0 projection)
                        (i/$ [0 1 2 3 4 5] 1 projection)
                        :title "hhv6 family"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [6 7 8 9 10 11 12 13]  0 projection)
                      (i/$ [6 7 8 9 10 11 12 13]  1 projection))
        (c/add-points (i/$ [14 15 16] 0 projection)
                      (i/$ [14 15 16] 1 projection))
        (c/add-points (i/$ [17 18] 0 projection)
                      (i/$ [17 18] 1 projection))
        (i/view))))

(defn HHV6-SVD-depth [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 1 2 9 10 14 15 16] 0 projection)
                        (i/$ [0 1 2 9 10 14 15 16] 1 projection)
                        :title "hhv6 depth"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [3 4 5 6 7 8 11 12 17 18]  0 projection)
                      (i/$ [3 4 5 6 7 8 11 12 17 18]  1 projection))
        (i/view))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Gneral description
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Gneral descriptive functions


(defn cutoff [file]
  (i/nrow (i/$where (i/$fn [minfr pi]
                           (and (> minfr 0.05)
                                (> pi 0.0))) file)))
(defn sum-cov [file]
  (i/sum (i/$ :depth file)))

(defn sum-pi [file]
  (i/sum (i/$ :pi (i/$where (i/$fn [minfr] (> minfr 0.05)) file))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Samples dataset EBC Taken from HCMV
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn sample-table [samples]
  (->> (i/dataset
        [ :sample            :player   :date]
        [[(samples :H37-Pa ) "Primary" 1]
         [(samples :H37-Pb ) "Primary" 1] ;;for graphs remove
         [(samples :H37-S2b) "Sibling" 2]
         [(samples :H37-S3 ) "Sibling" 3]
         [(samples :H37-Pc ) "Sibling" 1]
         [(samples :H37-S2a) "Sibling" 1]

         [(samples :H42-Ma ) "Mother"  1]
         [(samples :H42-Mb ) "Mother"  1]
         [(samples :H42-Pa ) "Primary" 1]
         [(samples :H42-Pb ) "Primary" 1]
         [(samples :H42-Pc ) "Primary" 1]
         [(samples :H42-S1a) "Sibling" 1]
         [(samples :H42-S1b) "Sibling" 1]
         [(samples :H42-S1c) "Sibling" 1]

         [(samples :H43-Pb ) "Primary" 1]
         [(samples :H43-S1a) "Sibling" 1]
         [(samples :H43-S1b) "Sibling" 1]

         [(samples :H72-M  ) "Mother"  1]
         [(samples :H72-Pa ) "Primary" 1]])
       
       (i/add-derived-column
        :name
        [:sample]
        #(first (i/$ :r_seq  %)))
       (i/add-derived-column
        :mean-cov
        [:sample]
        #(format "%.2f" (/ (sum-cov %)
                           (i/nrow  %))))
       (i/add-derived-column
        :cov>20
        [:sample]
        #(i/$where (i/$fn [cov] (> cov 20)) %))
       (i/add-derived-column
        :n-seg
        [:cov>20]
        #(if (> (i/nrow %) 0)
           (format "%.4f" (/ (double (cutoff %)) (i/nrow %))) 0.0))
       (i/add-derived-column
        :nuc-div
        [:cov>20]
        #(if (> (i/nrow %) 0)
           (format "%.4f" (/ (sum-pi %)  (i/nrow %))) 0.0))
       (i/add-derived-column
        :sfs
        [:cov>20]
        #(if (> (i/nrow %) 0) (p/bin-sfs 10 (da/get-synonymous %)) 0))
       (i/$ [:name :player :date :mean-cov :n-seg :nuc-div :sample :cov>20 :sfs])))
(def p-sample-table (memoize sample-table))


;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Simple visualizations

(defn show-one [sams [column title result]]
  (i/with-data sams
    (i/view (c/bar-chart :name column
                         :title title
                         :group-by :player
                         :x-label "Sample name"
                         :y-label result
                         :legend true
                         :vertical false))))

(defn show-box [sams [column title result]]
  (-> (c/box-plot column
                  :title title
                  :x-label "Sample name"
                  :y-label result
                  :legend true
                  :series-label "Primary"
                  :data (i/$where {:player {:$eq "Primary"}} sams))
      (c/add-box-plot column
                      :data (i/$where {:player {:$eq "Mother"}} sams)
                      :series-label "Mother")
      (c/add-box-plot column
                      :data (i/$where {:player {:$eq "Sibling"}} sams)
                      :series-label "Sibling")
            (i/view)))


(defn show-common [sams [column title result]]
  (i/with-data  (i/$rollup :mean column :player sams)
    (i/view (c/bar-chart :player column
                         :title title
                         :x-label "Sample name"
                         :y-label result
                         :legend true
                         :vertical false))))

(defn show-all [fnc]
  (let [sams (p-samples)]
    (map #(fnc sams %)

         [[:mean-cov "Mean coverage"                   "mean coverage"        ]
          [:nuc-div  "Nucleotide diversity"            "Nucleotide diversity" ] 
          [:n-seg    "Proportion of segregating sites" "Segregation propotion"]])))
#_(show-all show-one)
#_(show-all show-common)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;ALL SAMPLES descriptive ANALYSIS 

(defn run-all [fnc]
  " a function that accept funcinos returns a valur "
  (let [run_vec (i/$ :cov>20 (p-sample-table))]
    (double (/ (reduce + (map #(fnc %) run_vec))
               (reduce + (map #(i/nrow %) run_vec))))))

(defn stat-all []
  (println "mean coverage for all:           " (run-all sum-cov))
  (println "Nucleotide diversity:            " (run-all sum-pi ))
  (println "Segregating Sites per nucleotid: " (run-all cutoff)))
#_ (stat-all)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;SINGLE SAMPLE descriptive ANALYSIS

(defn summary-sfs [file]
  (let [synonymous (da/get-synonymous file)
        syn-sfs    (p/bin-sfs 10 synonymous)
        mean_cov   (/ (sum-cov file)      (i/nrow file))
        nuc_div    (/ (sum-pi file)       (i/nrow file))
        snp>nucs   (double (/ (cutoff file)       (i/nrow file)))]
    (println "\nSummary statistics for" (first (i/$ :r_seq file))":")
    (println "mean coverage:                    " mean_cov)
    (println "Nucleotide diversity:             " nuc_div)
    (println "Segregating Sites per nucleotid:  " snp>nucs)
    (println "Synonymous site frequency spectra:")
    (println (map first  syn-sfs))
    (println (map second syn-sfs))))

(defn stat-sample [fnc]
  (let [run_vec (i/$ :sample (p-samples))]
    (map #(fnc %) run_vec)))
#_(stat-sample summary-sfs)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET WINDOWED SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn win-100 [file]
  (p/m-slide-mean file :pi :pi_slide 100)) 

(defn les-win []
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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Non descriptive single sample
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Functions for understanding minor allele cutoff

(defn fr-dist [dep mfr file] 
  "tool for visualizing minor variant frequencies and FP + Depth"
  (i/$where (i/$fn [depth minfr]
                   (and (> depth dep)
                        (< minfr mfr)))  
            (i/$ [:ref-loc  :gfwd+ :gfwd- 
                  :CDS+     :CDS-  :ref :loc 
                  :depth :T :A  :C :G   :minfr :pi ] 
                 file)))

(defn all-freq-view [file1 file2 dep]
  (-> (c/xy-plot   :loc :minfr
                   :x-label "Position" :y-label "Minor variants frequency"
                   :title (str dep " minimal")
                   :data (fr-dist file1 dep 1.0)) 
      (c/add-lines :loc :minfr
                   :data (fr-dist file2 dep 1.0 ))
      (c/add-lines :loc  :minfr
                   :data (fr-dist file1 dep 0.1))
      (c/add-lines :loc  :minfr
                   :data (fr-dist file2 dep 0.1))
      (c/add-lines :loc  :minfr
                   :data (fr-dist file1 dep 0.003))
      (c/add-lines :loc  :minfr
                   :data (fr-dist file2 dep 0.003))
      (c/add-lines :loc 
                   (map #(/ % 100000) (i/$ :depth file1))
                   :data  file1)
      (c/add-lines :loc 
                   (map #(/ % 10000) (i/$ :depth file2))
                   :data  file2)
      (i/view)))   

(defn poisson-nonfilthered [dep mfr file]
  "get all poisson filtered data which is poistive under certain minor freq (mfr)"
  (i/$where (i/$fn [pi] (> pi 0.0)) (fr-dist file dep mfr)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Different Data slicing methods

(defn pi-chart [file]
  (->> file
       (i/$where (i/$fn [CDS+] (not= CDS+ "-")))
       (i/$ [:ref-loc :gfwd+ :gbwd+ :CDS+ :loc :pi
             :depth :maj+ :min+ :maj_aa+ :min_aa+ 
             :T     :A    :C    :G 
             :orf+ :majorf+ :minorf+])
       (i/$where (i/$fn [pi] (= pi 0.0)))))

(defn nonsyn-chart [file]
  (->> file
       (i/$where (i/$fn [CDS+] (not= CDS+ "-")))
       (i/$ [:ref-loc :gfwd+ :CDS+ :loc 
             :depth :maj+ :min+ :maj_aa+ :min_aa+ 
             :T     :A    :C    :G 
             :orf+ :majorf+ :minorf+])
       (i/$where (i/$fn [majorf+ minorf+] (not= majorf+ minorf+)))))

