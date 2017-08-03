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

