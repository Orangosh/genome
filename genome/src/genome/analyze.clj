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
           [genome.compare    :as gc ]
           [genome.view       :as v  ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(def home                    "/mnt/data/datafiles/")
(def input_file  (str home "concensus/merlin.gff3"))
(def output_file (str home "concensus/refset.inc" ))


(def L05-Pa  (str home"incanted_files/505-Pa.inc" ))
(def L05-M   (str home"incanted_files/505-M.inc"  ))

(def L19-Pb  (str home "incanted_files/519-Pb.inc"))
(def L19-Pc  (str home"incanted_files/519-Pc.inc" ))
(def L19-Pd  (str home"incanted_files/519-Pd.inc" ))
(def L19-S1a (str home"incanted_files/519-S1a.inc"))

(def L20-Pa  (str home"incanted_files/520-Pa.inc" ))
(def L20-Pb  (str home"incanted_files/520-Pb.inc" ))
(def L20-Pc  (str home"incanted_files/520-Pc.inc" ))
(def L20-S1  (str home"incanted_files/520-S1.inc" )) 
(def L20-S1a (str home"incanted_files/520-S1a.inc"))
  
(def L79-Pa  (str home"incanted_files/579-Pa.inc" ))
(def L79-Pb  (str home"incanted_files/579-Pb.inc" ))
(def L79-M   (str home"incanted_files/579-M.inc"  ))
(def L79-S1a (str home"incanted_files/579-S1a.inc"))
(def L79-S1b (str home"incanted_files/579-S1b.inc"))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-set [file cov]
  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [cov_p] (< cov cov_p)))))
  (def m-get-set (memoize get-set))

(defn les-sets []
  (def S05-Pa  (m-get-set L05-Pa  20))
  (def S05-M   (m-get-set L05-M   20))

  (def S19-Pb  (m-get-set L19-Pb  20))
  (def S19-Pc  (m-get-set L19-Pc  20))
  (def S19-Pd  (m-get-set L19-Pd  20))
  (def S19-S1a (m-get-set L19-S1a 20))

  (def S20-Pa  (m-get-set L20-Pa  20))
  (def S20-Pb  (m-get-set L20-Pb  20))
  (def S20-Pc  (m-get-set L20-Pc  20))
  (def S20-S1  (m-get-set L20-S1  20)) 
  (def S20-S1a (m-get-set L20-S1a 20))
  
  (def S79-Pa  (m-get-set L79-Pa  20))
  (def S79-Pb  (m-get-set L79-Pb  20))
  (def S79-M   (m-get-set L79-M   20))
  (def S79-S1a (m-get-set L79-S1a 20))
  (def S79-S1b (m-get-set L79-S1b 20)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET WINDOWED SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn win-100 [file]
  (p/m-slide-mean file :pi_pois :pi_slide 100)) 

(defn les-win []
  (def W05-Pa  (win-100 (m-get-set L05-Pa )))
  (def W05-M   (win-100 (m-get-set L05-M  )))

  (def W19-Pb  (win-100 (m-get-set L19-Pb )))
  (def W19-Pc  (win-100 (m-get-set L19-Pc )))
  (def W19-Pd  (win-100 (m-get-set L19-Pd )))
  (def W19-S1a (win-100 (m-get-set L19-S1a)))

  (def W20-Pa  (win-100 (m-get-set L20-Pa )))
  (def W20-Pb  (win-100 (m-get-set L20-Pb )))
  (def W20-Pc  (win-100 (m-get-set L20-Pc )))
  (def W20-S1  (win-100 (m-get-set L20-S1 )))
  (def W20-S1a (win-100 (m-get-set L20-S1a)))
  
  (def W79-Pa  (win-100 (m-get-set L79-Pa )))
  (def W79-Pb  (win-100 (m-get-set L79-Pb )))
  (def W79-M   (win-100 (m-get-set L79-M  )))
  (def W79-S1a (win-100 (m-get-set L79-S1a)))
  (def W79-S1b (win-100 (m-get-set L79-S1b))))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SINGLE SAMPLE ANALYSIS- DESCRIPTIVE STATISTICS  + SFS 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn sumary-sfs [file]
  (let [synonymous (da/get-synonymous file)
        syn-sfs (p/bin-sfs 10 synonymous)]
  (println "SUMMARY STATISTICS:")
  (gs/stat-report file)
  (println "Synonymous site frequency spectra")
  (println (map first  syn-sfs))
  (println (map second syn-sfs))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SINGLE SAMPLE ANALYSIS- GET COMMON FEATURES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn pi-chart [file]
  (->> file
       (i/$where (i/$fn [CDS+] (not= CDS+ "-")))
       (i/$ [:ref-loc :gene+ :CDS+ :loc :pi_pois
             :cov_p :maj_p+ :min_p+ :maj_aa+ :min_aa+ 
             :Tpois :Apois :Cpois :Gpois 
             :orf+ :majorf+ :minorf+])
       (i/$where (i/$fn [pi_pois] (= pi_pois 0.0)))))

(defn nonsyn-chart [file]
  (->> file
       (i/$where (i/$fn [CDS+] (not= CDS+ "-")))
       (i/$ [:ref-loc :gene+ :CDS+ :loc 
             :cov_p :maj_p+ :min_p+ :maj_aa+ :min_aa+ 
             :Tpois :Apois :Cpois :Gpois 
             :orf+ :majorf+ :minorf+])
       (i/$where (i/$fn [majorf+ minorf+] (not= majorf+ minorf+)))))


(import '(org.jfree.chart StandardChartTheme)
        '(org.jfree.chart.plot DefaultDrawingSupplier)
        '(java.awt Color))

(def all-red-theme 
  (doto (StandardChartTheme/createJFreeTheme)
    (.setDrawingSupplier
     (proxy [DefaultDrawingSupplier] []
       (getNextPaint [] Color/blue)))))

(defn get-gene [function col file]
  (let [general (->> file
                     (i/$ col)
                     frequencies
                     (map vec)
                     vec
                     (i/dataset [col :refsum]))
        specific (->> (function file)
                      (i/$ col)
                      frequencies
                      (map vec)
                      vec
                      (i/dataset [col :sum]))
        spec&gen (i/$join [col col] general specific)]
    (->> spec&gen
         (i/add-derived-column
          :ratio
          [:sum :refsum]
          #(double (/ %1 %2)))
         (i/$order :ratio :desc))))


(defn view-gene [function col filename file]
  (-> (c/bar-chart
       col
       :ratio
       :title filename 
       :legend true
       :data (i/$ (range 0 10) :all (get-gene file function col)))
      (c/set-theme all-red-theme)
      i/view))

(defn pi-view
  [file]
  (-> (c/xy-plot
       :loc
       :ratio
       :x-label "Position"
       :y-label "Ratio"
       :title "Comparde"
       :data file) 
      (c/add-lines
       :loc :rations
       :data file)
      (i/view)))   




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;COMPARISON ANALYSES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn seq-compare [comp_type file1 file2 ]
  (case comp_type
    "n-var" (gc/nuc-variants     file1 file2)
    "a-var" (gc/aa-variants      file1 file2)
    "all"   (gc/allele-change    file1 file2)
    "div"   (gc/diversity-change file1 file2)))
(def m-compare (memoize seq-compare))

(defn clean-compare [file1 file2 comp_type]
  "Adds only annotation rows that have a common :loc value with seq_dataset"
  (->> (m-compare file1 file2 comp_type)
        (i/$ [:loc  :gene+ :gene-  :mRNA+  :mRNA-
              :cov1 :cov2  :p-sus1 :p-sus2 :n-sus1 :n-sus2])
        (i/$where (i/$fn [p-sus1 p-sus2] (not= p-sus1 p-sus2)))
        (i/$where (i/$fn [cov1 cov2] (and (< 20 cov1) (< 20 cov2))))))
