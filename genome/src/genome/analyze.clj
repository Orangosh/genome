(ns genome.analyze
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [genome.consvar    :as cv ]
           [genome.pop        :as p  ]
           [genome.compare    :as gc ]
           [genome.view       :as v  ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def input_file       "/home/yosh/datafiles/genes/merlin.gff3")
(def output_file      "/home/yosh/datafiles/genes/merlin.inc" )


(def L05-Pa  "/home/yosh/datafiles/incanted_files/505-Pa.inc" )
(def L05-M   "/home/yosh/datafiles/incanted_files/505-M.inc"  )

(def L19-Pb  "/home/yosh/datafiles/incanted_files/519-Pb.inc" )
(def L19-Pc  "/home/yosh/datafiles/incanted_files/519-Pc.inc" )
(def L19-Pd  "/home/yosh/datafiles/incanted_files/519-Pd.inc" )
(def L19-S1a "/home/yosh/datafiles/incanted_files/519-S1a.inc")

(def L20-Pa  "/home/yosh/datafiles/incanted_files/520-Pa.inc" )
(def L20-Pb  "/home/yosh/datafiles/incanted_files/520-Pb.inc" )
(def L20-Pc  "/home/yosh/datafiles/incanted_files/520-Pc.inc" )
(def L20-S1  "/home/yosh/datafiles/incanted_files/520-S1.inc" ) 
(def L20-S1a "/home/yosh/datafiles/incanted_files/520-S1a.inc")
  
(def L79-Pa  "/home/yosh/datafiles/incanted_files/579-Pa.inc" )
(def L79-Pb  "/home/yosh/datafiles/incanted_files/579-Pb.inc" )
(def L79-M   "/home/yosh/datafiles/incanted_files/579-M.inc"  )
(def L79-S1a "/home/yosh/datafiles/incanted_files/579-S1a.inc")
(def L79-S1b "/home/yosh/datafiles/incanted_files/579-S1b.inc")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-set [file]
  "open an csv.inc file"
  (ii/read-dataset file :header true))
(def m-get-set (memoize get-set))

(def S05-Pa  (m-get-set L05-Pa ))
(def S05-M   (m-get-set L05-M  ))

(def S19-Pb  (m-get-set L19-Pb ))
(def S19-Pc  (m-get-set L19-Pc ))
(def S19-Pd  (m-get-set L19-Pd ))
(def S19-S1a (m-get-set L19-S1a))

(def S20-Pa  (m-get-set L20-Pa ))
(def S20-Pb  (m-get-set L20-Pb ))
(def S20-Pc  (m-get-set L20-Pc ))
(def S20-S1  (m-get-set L20-S1 )) 
(def S20-S1a (m-get-set L20-S1a))
  
(def S79-Pa  (m-get-set L79-Pa ))
(def S79-Pb  (m-get-set L79-Pb ))
(def S79-M   (m-get-set L79-M  ))
(def S79-S1a (m-get-set L79-S1a))
(def S79-S1b (m-get-set L79-S1b))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET WINDOWED SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn win-1000 [file]
  (p/m-slide-mean file :pi_pois :pi_slide 1000)) 

(def W05-Pa  (win-1000 (m-get-set L05-Pa )))
(def W05-M   (win-1000 (m-get-set L05-M  )))

(def W19-Pb  (win-1000 (m-get-set L19-Pb )))
(def W19-Pc  (win-1000 (m-get-set L19-Pc )))
(def W19-Pd  (win-1000 (m-get-set L19-Pd )))
(def W19-S1a (win-1000 (m-get-set L19-S1a)))

(def W20-Pa  (win-1000 (m-get-set L20-Pa )))
(def W20-Pb  (win-1000 (m-get-set L20-Pb )))
(def W20-Pc  (win-1000 (m-get-set L20-Pc )))
(def W20-S1  (win-1000 (m-get-set L20-S1 )))
(def W20-S1a (win-1000 (m-get-set L20-S1a)))
  
(def W79-Pa  (win-1000 (m-get-set L79-Pa )))
(def W79-Pb  (win-1000 (m-get-set L79-Pb )))
(def W79-M   (win-1000 (m-get-set L79-M  )))
(def W79-S1a (win-1000 (m-get-set L79-S1a)))
(def W79-S1b (win-1000 (m-get-set L79-S1b)))

(v/look :pi_slide S79-Pb1000 S79-M1000 S79-S1a1000)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SINGLE SAMPLE ANALYSIS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn sumary-sfs [file]
  (let [synonymous (da/get-synonymous file)
        syn-sfs (p/bin-sfs 10 synonymous)]
  (println "SUMMARY STATISTICS:")
  (gs/stat-report file)
  (println "Synonymous site frequency spectra")
  (println (map first  syn-sfs))
  (println (map second syn-sfs))

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;COMPARISON ANALYSES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn seq-compare [file1 file2 comp_type]
  (case comp_type
    "n-var" (gc/nuc-variants   file1 file2)
    "a-var" (gc/aa-variants    file1 file2)
    "all" (gc/allele-change    file1 file2)
    "div" (gc/diversity-change file1 file2)))
(def m-compare (memoize seq-compare))

(defn clean-compare [file1 file2 comp_type]
  "Adds only annotation rows that have a common :loc value with seq_dataset"
  (->> (m-compare file1 file2 comp_type)
        (i/$ [:loc  :gene+ :gene-  :mRNA+  :mRNA-
              :cov1 :cov2  :p-sus1 :p-sus2 :n-sus1 :n-sus2])
        (i/$where (i/$fn [p-sus1 p-sus2] (not= p-sus1 p-sus2)))
        (i/$where (i/$fn [cov1 cov2] (and (< 20 cov1) (< 20 cov2))))))
