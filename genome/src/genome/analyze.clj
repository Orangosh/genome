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
           [clojure.data.csv  :as csv]
           [genome.view       :as v  ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;
;; For PC

;;(def home                  "/home/yosh/datafiles/")
;;(def input_file  (str home "genes/merlin.gff3"   ))
;;(def output_file (str home "genes/merlin.inc"    ))

;;;;;;;;;;;;;;;;;;;;;;;
;; For Server

(def home                  "/mnt/data/datafiles/"  )
(def input_file  (str home "concensus/merlin.gff3"))
(def output_file (str home "concensus/refset.inc" ))


 
(def L05-Pa  (str home "incanted_files/505-Pa.inc" ))
(def L05-M   (str home "incanted_files/505-M.inc"  ))

(def L19-Pb  (str home "incanted_files/519-Pb.inc" ))
(def L19-Pc  (str home "incanted_files/519-Pc.inc" ))
(def L19-Pd  (str home "incanted_files/519-Pd.inc" ))
(def L19-S1a (str home "incanted_files/519-S1a.inc"))

(def L20-Pa  (str home "incanted_files/520-Pa.inc" ))
(def L20-Pb  (str home "incanted_files/520-Pb.inc" ))
(def L20-Pc  (str home "incanted_files/520-Pc.inc" ))
(def L20-S1  (str home "incanted_files/520-S1.inc" )) 
(def L20-S1a (str home "incanted_files/520-S1a.inc"))
  
(def L79-Pa  (str home "incanted_files/579-Pa.inc" ))
(def L79-Pb  (str home "incanted_files/579-Pb.inc" ))
(def L79-M   (str home "incanted_files/579-M.inc"  ))
(def L79-S1a (str home "incanted_files/579-S1a.inc"))
(def L79-S1b (str home "incanted_files/579-S1b.inc"))


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
(def S20-S1  (m-get-set L20-S1  0)) 
(def S20-S1a (m-get-set L20-S1a 0))
  
(def S79-Pa  (m-get-set L79-Pa  0))
(def S79-Pb  (m-get-set L79-Pb  0))
(def S79-M   (m-get-set L79-M   0))
(def S79-S1a (m-get-set L79-S1a 0))
(def S79-S1b (m-get-set L79-S1b 0))

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
  (def W20-S1  (win-100 (m-get-set L20-S1  20)))
  (def W20-S1a (win-100 (m-get-set L20-S1a 20)))
  
  (def W79-Pa  (win-100 (m-get-set L79-Pa  20)))
  (def W79-Pb  (win-100 (m-get-set L79-Pb  20)))
  (def W79-M   (win-100 (m-get-set L79-M   20)))
  (def W79-S1a (win-100 (m-get-set L79-S1a 20)))
  (def W79-S1b (win-100 (m-get-set L79-S1b 20))))




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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Functions for understanding minor allele cutoff

(defn fr-dist [dep mfr file] 
  "toll for visualizing minor variant frequencies and FP + Depth"
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
       (i/$ [:ref-loc :gfwd+ :CDS+ :loc :pi
             :depth :maj+ :min+ :maj_aa+ :min_aa+ 
             :T     :A    :C    :G 
             :orf+ :majorf+ :minorf+])
       (i/$where (i/$fn [pi_pois] (= pi_pois 0.0)))))

(defn nonsyn-chart [file]
  (->> file
       (i/$where (i/$fn [CDS+] (not= CDS+ "-")))
       (i/$ [:ref-loc :gfwd+ :CDS+ :loc 
             :depth :maj+ :min+ :maj_aa+ :min_aa+ 
             :T     :A    :C    :G 
             :orf+ :majorf+ :minorf+])
       (i/$where (i/$fn [majorf+ minorf+] (not= majorf+ minorf+)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Get genes

(import '(org.jfree.chart StandardChartTheme)
        '(org.jfree.chart.plot DefaultDrawingSupplier)
        '(java.awt Color))

(def all-red-theme 
  (doto (StandardChartTheme/createJFreeTheme)
    (.setDrawingSupplier
     (proxy [DefaultDrawingSupplier] []
       (getNextPaint [] Color/blue)))))

(defn get-gene [function col file]
  "function is any filter/ col is usualy :gfwd+/-"
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

(defn filtre [file] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file
       (i/$where (i/$fn [depth1 minfr1 depth2 minfr2
                         depth3 minfr3 depth4 minfr4]
                        (and (not= nil depth1) (not= nil depth2)
                             (not= nil depth3) (not= nil depth4))))
       (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2
                         majorf+3 minorf+3 majorf+4 minorf+4]
                        (and (not= majorf+1 minorf+1)
                             (not= majorf+2 minorf+2)
                             (not= majorf+3 minorf+3)
                             (not= majorf+4 minorf+4))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3 minfr4]
                        (and (> minfr1 0.003) (> minfr2 0.003)
                             (> minfr3 0.003) (> minfr4 0.003))))
       (i/$where (i/$fn [depth1 depth2 depth3 depth4]
                        (and (> depth1 20.0 ) (> depth2 20.0 )
                             (> depth3 20.0 ) (> depth4 20.0 ))))
       (i/$where (i/$fn [CDS+ CDS-]
                        (or  (not= CDS+ "-" ) )))
       (i/$ [:ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS-
             :minfr1 :minfr2 :minfr3 :minfr4])))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PREPARING DATA FOR d3/CIRCOS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;
;;Anotations

(defn get-locations [col postrand? begining? file]
  "gets a set with all begining or end of a positive or netetive column"
  (let [gene     (if postrand? :gfwd+ :gfwd-)
        name-vec (->> file
                      (i/$where {col {:$ne "-"}})
                      (i/$ gene)
                      distinct)]
    (if begining?     
      (map #(i/$ (range 0 1) :all
                 (i/$where {col {:$ne "-"}}
                           (i/$where {gene {:$eq %}} file)))
           name-vec)
      (map #(i/$ (range 0 1) :all
                 (i/$order
                  :loc :desc
                  (i/$where {col {:$ne "-"}}
                            (i/$where {gene {:$eq %}} file))))
           name-vec))))

(defn get-genes [col postrand? begining? file ]
  (->>(get-locations col postrand? begining? file)
      (apply i/conj-rows)))
  (def m-get-genes (memoize get-genes))

(defn get-range-set [col postrand? file]
  (let [file+ (m-get-genes col postrand? true  file)
        file- (m-get-genes col postrand? false file)
        ncol (i/$ :loc file-)]
    (->> file+
         (i/add-column
          :end
          ncol)
         (i/rename-cols {:loc :start
                         :gfwd+ :gfwdp
                         :gfwd- :gfwdn
                         :CDS+  :CDSp
                         :CDS-  :CDSn})
         (i/$ [:r_seq  col   :start :end
               :gfwdp :gfwdn :CDSp  :CDSn]))))
(def m-get-range-set (memoize get-range-set))


;;;;;;;;;;;;;;;;;;
;;Diversity

(defn circos-database [file]
  (->> file
       (i/rename-cols {:r_seq   :sample
                       :loc     :position1
                       :pi_slide :value})
       (i/add-derived-column
        :position2
        [:position1]
        #(+ % 1))
       (i/$ [:sample :position1 :position2 :value])))

(defn circosing1 []
  (let [ An19-Pb (m-get-set L19-Pb 0)] 
    [[(m-get-range-set :gfwd+ true  An19-Pb)
      "/home/yosh/Software/git/visual/circos/resources/public/data/genepos.csv"]
     [(m-get-range-set :gfwd- false An19-Pb)
      "/home/yosh/Software/git/visual/circos/resources/public/data/geneneg.csv"]
     [(m-get-range-set :CDS+  true  An19-Pb)
      "/home/yosh/Software/git/visual/circos/resources/public/data/CDSpos.csv"]
     [(m-get-range-set :CDS-  false An19-Pb)
      "/home/yosh/Software/git/visual/circos/resources/public/data/CDSneg.csv"]]))

(defn circosing2 []     
  [[W05-Pa
    "/home/yosh/Software/git/visual/circos/resources/public/data/W05-Pa.csv" ]
   [W05-M
    "/home/yosh/Software/git/visual/circos/resources/public/data/W05-M.csv"  ]

   [W19-Pb
    "/home/yosh/Software/git/visual/circos/resources/public/data/W19-Pb.csv" ]
   [W19-Pc
    "/home/yosh/Software/git/visual/circos/resources/public/data/W19-Pc.csv" ]
   [W19-Pd
    "/home/yosh/Software/git/visual/circos/resources/public/data/W19-Pd.csv" ]
   [W19-S1a
    "/home/yosh/Software/git/visual/circos/resources/public/data/W19-S1a.csv"]])
(defn circosing3 []          
  [[W20-Pa
    "/home/yosh/Software/git/visual/circos/resources/public/data/W20-Pa.csv" ]
   [W20-Pb
    "/home/yosh/Software/git/visual/circos/resources/public/data/W20-Pb.csv" ]
   [W20-Pc
    "/home/yosh/Software/git/visual/circos/resources/public/data/W20-Pc.csv" ]
   [W20-S1
    "/home/yosh/Software/git/visual/circos/resources/public/data/W20-S1.csv" ]
   [S19-S1a
    "/home/yosh/Software/git/visual/circos/resources/public/data/W20-S1a.csv"]
   [W79-Pa
    "/home/yosh/Software/git/visual/circos/resources/public/data/W79-Pa.csv" ]
   [W79-Pb
    "/home/yosh/Software/git/visual/circos/resources/public/data/W79-Pb.csv" ]
   [W79-M
    "/home/yosh/Software/git/visual/circos/resources/public/data/W79-M.csv"  ]
   [W79-S1a
    "/home/yosh/Software/git/visual/circos/resources/public/data/W79-S1a.csv"]
   [W79-S1b
    "/home/yosh/Software/git/visual/circos/resources/public/data/W79-S1b.csv"]])


(defn circos [circosing]
  (let [[file_in file_out] circosing
        file_inc           (circos-database file_in)]
    (with-open [f-out (io/writer file_out)]
      (csv/write-csv f-out [(map name (i/col-names file_inc))])
      (csv/write-csv f-out (i/to-list file_inc)))))

(defn cir-ann [circosing]
  (let [[file_in file_out] circosing]
    (with-open [f-out (io/writer file_out)]
      (csv/write-csv f-out [(map name (i/col-names file_in))])
      (csv/write-csv f-out (i/to-list file_in)))))

(defn map-circos-sets [circosing]
  (map #(circos %) circosing))
