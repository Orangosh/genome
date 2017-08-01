(ns genome.analyze
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
;;(def home "/home/yosh/datafiles/incanted_files/")

;;;;;;;;;;;;;;;;;;;;;;;
;; For Server
(def home "/mnt/data/datafiles/incanted_files/")


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
(def L79-S1b (str home "579-S1b.inc"))

(defn bm-loc []
  (def L1  (str home "S1.inc"  ))
  (def L10 (str home "S10.inc" ))
  (def L11 (str home "S11.inc" ))
  (def L12 (str home "S12.inc" ))
  (def L13 (str home "S13.inc" ))
  (def L14 (str home "S14.inc" ))
  (def L15 (str home "S15.inc" ))
  (def L16 (str home "S16.inc" ))
  (def L17 (str home "S17.inc" ))
  (def L18 (str home "S18.inc" ))
  (def L23 (str home "S23.inc" ))
  (def L24 (str home "S24.inc" ))
  (def L25 (str home "S25.inc" ))
  (def L26 (str home "S26.inc" ))
  (def L27 (str home "S27.inc" ))
  (def L28 (str home "S28.inc" ))
  (def L29 (str home "S29.inc" ))
  (def L30 (str home "S30.inc" ))
  (def L9  (str home "S9.inc"  )))


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

(defn bm-sets []
  (def S1  (m-get-set L1  0))
  (def S10 (m-get-set L10 0))

  (def S11 (m-get-set L11 0))
  (def S12 (m-get-set L12 0))
  (def S13 (m-get-set L13 0))
  (def S14 (m-get-set L14 0))

  (def S15 (m-get-set L15 0))
  (def S16 (m-get-set L16 0))
  (def S17 (m-get-set L17 0))
  (def S18 (m-get-set L18 0)) 
  (def S23 (m-get-set L23 0))
  
  (def S24 (m-get-set L24 0))
  (def S25 (m-get-set L25 0))
  (def S26 (m-get-set L26 0))
  (def S27 (m-get-set L27 0))
  (def S28 (m-get-set L28 0))
  (def S29 (m-get-set L29 0))
  (def S30 (m-get-set L30 0))
  (def S9  (m-get-set L9  0)))


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
        [[S05-Pa  "Primary" 1]
         [S05-M   "Mother"  1];;for graphs remove
         [S19-Pb  "Primary" 1]
         [S19-Pc  "Primary" 2]
         [S19-Pd  "Primary" 3]
         [S19-S1a "Sibling" 1]
         [S20-Pa  "Primary" 1]
         [S20-Pb  "Primary" 2]
         [S20-Pc  "Primary" 3]
         [S20-S1a "Sibling" 1]
         [S20-S1b "Sibling" 2]
         [S79-Pa  "Primary" 1]
         [S79-Pb  "Primary" 2]
         [S79-M   "Mother"  0]
         [S79-S1a "Sibling" 1]
         [S79-S1b "Sibling" 2]])
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
  (let [run_vec (i/$ :cov>20 (p-samples))]
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

