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

(defn variants-comparisson [file1 file2]
  (let [seq1 (m-get-set file1)
        seq2 (m-get-set file2)]
    (->> (gc/variants seq1 seq2))))
  (def m-variants-comparisson (memoize variants-comparisson))

(defn annotate-compare [seq_dataset annotation]
  "Adds only annotation rows that have a common :loc value with seq_dataset"
  (->>  (i/$join [:loc :loc] seq_dataset annotation)
        (i/$ [:loc :gene+ :gene- :mRNA+ :mRNA-   :cov1 :cov2
              :p-sus1     :p-sus2       :n-sus1  :n-sus2
              :aa_fwd1    :aa_fwd2      :aa_rev1 :aa_rev2])
        (i/$where (i/$fn [p-sus1 p-sus2] (not= p-sus1 p-sus2)))
        (i/$where (i/$fn [cov1 cov2] (and (< 20 cov1) (< 20 cov2))))))
(def m-annotate-compare (memoize annotate-compare))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def input_file  "/home/yosh/datafiles/genes/merlin.gff3")
(def output_file "/home/yosh/datafiles/genes/merlin.inc" )

(def S79-Pa "/home/yosh/datafiles/incanted_files/579-Pa.inc" )
(def S79-Pb "/home/yosh/datafiles/incanted_files/579-Pb.inc" )
(def S79-M  "/home/yosh/datafiles/incanted_files/579-M.inc"  )
(def S19-Pb "/home/yosh/datafiles/incanted_files/519-Pb.inc" )
(def S19-Pc "/home/yosh/datafiles/incanted_files/519-Pc.inc" )
(def S19-Pd "/home/yosh/datafiles/incanted_files/519-Pd.inc" )
(def S19-S  "/home/yosh/datafiles/incanted_files/519-S1a.inc")

(def anncompared (m-annotate-compare (m-get-annotation input_file)
                                   (m-variants-comparisson S19-Pc S19-Pb)))