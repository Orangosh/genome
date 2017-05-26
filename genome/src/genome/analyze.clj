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
(def home                  "/home/yosh/datafiles/")
(def input_file  (str home "genes/merlin.gff3"  ))
(def output_file (str home "genes/merlin.inc"   ))


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
       (i/$where (i/$fn [cov_p] (< cov cov_p)))))
  (def m-get-set (memoize get-set))

(defn les-sets [])
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
(def S79-S1b (m-get-set L79-S1b 20))

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

(defn clean-proc [comp_file]
  "Adds only annotation rows that have a common :loc value with seq_dataset"
  (->> comp_file
       (i/$where (i/$fn [CDS+] (not= CDS+ "-")))
       (i/$where (i/$fn [majorf+1 majorf+2 minorf+1 minorf+2]
                        (and (not= majorf+2 minorf+2)
                             (=    minorf+2 majorf+1)
                             (not= majorf+1 "-")
                             (not= majorf+2 "-")
                             (not= minorf+1 "-")
                             (not= minorf+2 "-"))))
       (i/$ [:ref-loc1 :loc :gene+
             :cov1 :A1 :T1 :G1 :C1
             :cov2 :A2 :T2 :G2 :C2
             :majorf+1 :majorf+2 :minorf+1 :minorf+2])))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PREPARING DATA FOR d3/CIRCOS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;
;;Anotations

(defn get-locations [col postrand? begining? file]
  "gets a set with all begining or end of a positive or netetive column"
  (let [gene     (if postrand? :gene+ :gene-)
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
         (i/rename-cols {:loc :start})
         (i/$ [:r_seq :start :end
               :gene+ :gene- col]))))

;;;;;;;;;;;;;;;;;;
;;Diversity

(defn circos-database [file]
  (->> file
       (i/$ [:r_seq :loc :pi_pois])
       (i/rename-cols {:r_seq   :sample
                       :loc     :position1
                       :pi_pois :value})
       (i/add-derived-column
        :position2
        [:position1]
        (+ % 1))))

(def circosing
  [S19-Pb "/home/yosh/Software/git/circos/resources/public/data/S19-Pb.csv"])

(defn circos [circosing]
  (let [[file_in file_out] circosing
        file_inc           (circos-database file_in)]
    (with-open [f-out (io/writer file_out)]
      (csv/write-csv f-out [(map name (i/col-names file_inc))])
      (csv/write-csv f-out (i/to-list file_inc)))))

