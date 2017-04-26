(ns genome.analyze
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [genome.consvar    :as cv ]
           [genome.compare    :as gc ]
           [genome.annotate   :as ga ]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PROCESSES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-annotation [input_file]
  "creates an annotation"
  (->> (ga/gff3>dataset input_file)
       (i/$ [:loc :gene+ :gene- :mRNA+ :mRNA-]) ))
  (def m-get-annotation (memoize get-annotation))

(defn get-set [file]
  "open an csv.inc file"
  (ii/read-dataset file :header true))
  (def m-get-set (memoize get-set))

(defn annotate-compare [file1 file2 comp_type input_file]
  (let [seq1 (m-get-set file1)
        seq2 (m-get-set file2)
        annotation (m-get-annotation input_file)
        typed_set (case comp_type
                    "var" (gc/variants         seq1 seq2)
                    "all" (gc/allele-change    seq1 seq2)
                    "div" (gc/diversity-change seq1 seq2))]
    (i/$join [:loc :loc] annotation typed_set )))
(def m-annotate-compare (memoize annotate-compare))

(defn clean-compare [file1 file2 comp_type annotation]
  "Adds only annotation rows that have a common :loc value with seq_dataset"
  (->> (m-annotate-compare file1 file2 comp_type annotation)
        (i/$ [:loc  :gene+ :gene-  :mRNA+  :mRNA-
              :cov1 :cov2  :p-sus1 :p-sus2 :n-sus1 :n-sus2])
        (i/$where (i/$fn [p-sus1 p-sus2] (not= p-sus1 p-sus2)))
        (i/$where (i/$fn [cov1 cov2] (and (< 20 cov1) (< 20 cov2))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def input_file  "/home/yosh/datafiles/genes/merlin.gff3")
(def output_file "/home/yosh/datafiles/genes/merlin.inc" )


(def S05-Pa  "/home/yosh/datafiles/incanted_files/505-Pa.inc" )
(def S05-M   "/home/yosh/datafiles/incanted_files/505-M.inc"  )

(def S19-Pb  "/home/yosh/datafiles/incanted_files/519-Pb.inc" )
(def S19-Pc  "/home/yosh/datafiles/incanted_files/519-Pc.inc" )
(def S19-Pd  "/home/yosh/datafiles/incanted_files/519-Pd.inc" )
(def S19-S1a "/home/yosh/datafiles/incanted_files/519-S1a.inc")

(def S20-Pa  "/home/yosh/datafiles/incanted_files/520-Pa.inc" )
(def S20-Pb  "/home/yosh/datafiles/incanted_files/520-Pb.inc" )
(def S20-Pc  "/home/yosh/datafiles/incanted_files/520-Pc.inc" )
(def S20-S1  "/home/yosh/datafiles/incanted_files/520-S1.inc" ) 
(def S20-S1a "/home/yosh/datafiles/incanted_files/520-S1a.inc")
  
(def S79-Pa  "/home/yosh/datafiles/incanted_files/579-Pa.inc" )
(def S79-Pb  "/home/yosh/datafiles/incanted_files/579-Pb.inc" )
(def S79-M   "/home/yosh/datafiles/incanted_files/579-M.inc"  )
(def S79-S1a "/home/yosh/datafiles/incanted_files/579-S1a.inc")
(def S79-S1b "/home/yosh/datafiles/incanted_files/579-S1b.inc")


