(ns genome.getseqs
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


(defn get-loc []
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

  (def L9      (str home "bmseq/incanted_files/S9.inc"     ))

  ;;cmv-sets
  (def L05-Pa  (str home "bmseq/incanted_files/505-Pa.inc" ))
  (def L05-M   (str home "bmseq/incanted_files/505-M.inc"  ))

  (def L19-Pb  (str home "bmseq/incanted_files/519-Pb.inc" ))
  (def L19-Pc  (str home "bmseq/incanted_files/519-Pc.inc" ))
  (def L19-Pd  (str home "bmseq/incanted_files/519-Pd.inc" ))
  (def L19-S1a (str home "bmseq/incanted_files/519-S1a.inc"))

  (def L20-Pa  (str home "bmseq/incanted_files/520-Pa.inc" ))
  (def L20-Pb  (str home "bmseq/incanted_files/520-Pb.inc" ))
  (def L20-Pc  (str home "bmseq/incanted_files/520-Pc.inc" ))
  (def L20-S1  (str home "bmseq/incanted_files/520-S1.inc" ))
  (def L20-S1a (str home "bmseq/incanted_files/520-S1a.inc"))

  (def L79-Pa  (str home "bmseq/incanted_files/579-Pa.inc" ))
  (def L79-Pb  (str home "bmseq/incanted_files/579-Pb.inc" ))
  (def L79-M   (str home "bmseq/incanted_files/579-M.inc"  ))
  (def L79-S1a (str home "bmseq/incanted_files/579-S1a.inc"))
  (def L79-S1b (str home "bmseq/incanted_files/579-S1b.inc"))

  ;;hhv6-sets 
  (def L37-Pa  (str home "hhv6/incanted_files/537-Pa.inc"  ))
  (def L37-S2b (str home "hhv6/incanted_files/537-S2b.inc" ))        
  (def L37-Pb  (str home "hhv6/incanted_files/537-Pb.inc"  ))
  (def L37-S3  (str home "hhv6/incanted_files/537-S3.inc"  ))       
  (def L37-Pc  (str home "hhv6/incanted_files/537-Pc.inc"  ))            
  (def L37-S2a (str home "hhv6/incanted_files/537-S2a.inc" ))

  (def L42-Mb  (str home "hhv6/incanted_files/542-Mb.inc"  ))
  (def L42-Ma  (str home "hhv6/incanted_files/542-Ma.inc"  ))
  (def L42-Pa  (str home "hhv6/incanted_files/542-Pa.inc"  ))
  (def L42-Pb  (str home "hhv6/incanted_files/542-Pb.inc"  ))
  (def L42-Pc  (str home "hhv6/incanted_files/542-Pc.inc"  ))
  (def L42-S1a (str home "hhv6/incanted_files/542-S1a.inc" ))
  (def L42-S1b (str home "hhv6/incanted_files/542-S1b.inc" ))
  (def L42-S1c (str home "hhv6/incanted_files/542-S1c.inc" ))

  (def L43-Pb  (str home "hhv6/incanted_files/543-Pb.inc"  ))
  (def L43-S1a (str home "hhv6/incanted_files/543-S1a.inc" ))  
  (def L43-S1b (str home "hhv6/incanted_files/543-S1b.inc" ))

  (def L72-M   (str home "hhv6/incanted_files/572-M.inc"   ))
  (def L72-Pa  (str home "hhv6/incanted_files/572-Pa.inc"  ))

  ;; ebv-sets []
  (def L44-Pb  (str home "ebv/incanted_files/344-Pb.inc"   ))   
  (def L44-S2a (str home "ebv/incanted_files/344-S2.inc"   ))  
  (def L44-S2b (str home "ebv/incanted_files/344-S2b.inc"  ))  
  (def L44-S2c (str home "ebv/incanted_files/344-S2c.inc"  ))   
  (def L44-S3  (str home "ebv/incanted_files/344-S3.inc"   ))
  
  (def L25-Ma  (str home "ebv/incanted_files/525-Ma.inc"   ))  
  (def L25-Mc  (str home "ebv/incanted_files/525-Mc.inc"   ))  
  (def L25-Pa  (str home "ebv/incanted_files/525-Pa.inc"   ))  
  (def L25-S1a (str home "ebv/incanted_files/525-S1a.inc"  ))  
  (def L25-S1b (str home "ebv/incanted_files/525-S1b.inc"  ))  
  (def L25-S1c (str home "ebv/incanted_files/525-S1c.inc"  ))  
  (def L25-S1d (str home "ebv/incanted_files/525-S1d.inc"  ))  
  (def L25-S2a (str home "ebv/incanted_files/525-S2a.inc"  ))  
  (def L25-S2b (str home "ebv/incanted_files/525-S2b.inc"  ))  
  (def L25-S3a (str home "ebv/incanted_files/525-S3a.inc"  ))  
  (def L25-S3b (str home "ebv/incanted_files/525-S3b.inc"  ))  
  (def L25-S3c (str home "ebv/incanted_files/525-S3c.inc"  ))  

  (def L38-Ma  (str home "ebv/incanted_files/538-Ma.inc"   ))   
  (def L38-Pb  (str home "ebv/incanted_files/538-Pb.inc"   ))   
  (def L38-S1a (str home "ebv/incanted_files/538-S1a.inc"  ))  

  (def L40-Ma  (str home "ebv/incanted_files/540-Ma.inc"   ))  
  (def L40-Mb  (str home "ebv/incanted_files/540-Mb.inc"   ))  
  (def L40-Mc  (str home "ebv/incanted_files/540-Mc.inc"   ))  
  (def L40-Pa  (str home "ebv/incanted_files/540-Pa.inc"   ))
  (def L40-Pc  (str home "ebv/incanted_files/540-Pc.inc"   ))
  (def L40-S1a (str home "ebv/incanted_files/540-S1a.inc"  ))
  (def L40-S1b (str home "ebv/incanted_files/540-S1b.inc"  )))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-set [file cov]
  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))

(defn hcmv-sets []
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
  (def S79-S1b (m-get-set L79-S1b 0)))

(defn bm-sets []
  (def S1      (m-get-set L1      0))
  (def S10     (m-get-set L10     0))

  (def S11     (m-get-set L11     0))
  (def S12     (m-get-set L12     0))
  (def S13     (m-get-set L13     0))
  (def S14     (m-get-set L14     0))

  (def S15     (m-get-set L15     0))
  (def S16     (m-get-set L16     0))
  (def S17     (m-get-set L17     0))
  (def S18     (m-get-set L18     0))
  (def S23     (m-get-set L23     0))

  (def S24     (m-get-set L24     0))
  (def S25     (m-get-set L25     0))
  (def S26     (m-get-set L26     0))
  (def S27     (m-get-set L27     0))
  (def S28     (m-get-set L28     0))
  (def S29     (m-get-set L29     0))
  (def S30     (m-get-set L30     0))
  (def S9      (m-get-set L9      0)))

(defn hhv6-sets []
  (def H37-Pa  (m-get-set L37-Pa  0))
  (def H37-S2b (m-get-set L37-S2b 0))        
  (def H37-Pb  (m-get-set L37-Pb  0))
  (def H37-S3  (m-get-set L37-S3  0))       
  (def H37-Pc  (m-get-set L37-Pc  0))            
  (def H37-S2a (m-get-set L37-S2a 0))

  (def H42-Mb  (m-get-set L42-Mb  0))
  (def H42-Ma  (m-get-set L42-Ma  0))
  (def H42-Pa  (m-get-set L42-Pa  0))
  (def H42-Pb  (m-get-set L42-Pb  0))
  (def H42-Pc  (m-get-set L42-Pc  0))
  (def H42-S1a (m-get-set L42-S1a 0))
  (def H42-S1b (m-get-set L42-S1b 0))
  (def H42-S1c (m-get-set L42-S1c 0))

  (def H43-Pb  (m-get-set L43-Pb  0))
  (def H43-S1a (m-get-set L43-S1a 0))  
  (def H43-S1b (m-get-set L43-S1b 0))

  (def H72-M   (m-get-set L72-M   0))
  (def H72-Pa  (m-get-set L72-Pa  0)))

(defn ebv-sets []
  (def E44-Pb  (m-get-set L44-Pb  0))   
  (def E44-S2a (m-get-set L44-S2a 0))  
  (def E44-S2b (m-get-set L44-S2b 0))  
  (def E44-S2c (m-get-set L44-S2c 0))   
  (def E44-S3  (m-get-set L44-S3  0))
  
  (def E25-Ma  (m-get-set L25-Ma  0))  
  (def E25-Mc  (m-get-set L25-Mc  0))  
  (def E25-Pa  (m-get-set L25-Pa  0))  
  (def E25-S1a (m-get-set L25-S1a 0))  
  (def E25-S1b (m-get-set L25-S1b 0))  
  (def E25-S1c (m-get-set L25-S1c 0))  
  (def E25-S1d (m-get-set L25-S1d 0))  
  (def E25-S2a (m-get-set L25-S2a 0))  
  (def E25-S2b (m-get-set L25-S2b 0))  
  (def E25-S3a (m-get-set L25-S3a 0))  
  (def E25-S3b (m-get-set L25-S3b 0))  
  (def E25-S3c (m-get-set L25-S3c 0))  

  (def E38-Ma  (m-get-set L38-Ma  0))   
  (def E38-Pb  (m-get-set L38-Pb  0))   
  (def E38-S1a (m-get-set L38-S1a 0))  

  (def E40-Ma  (m-get-set L40-Ma  0))  
  (def E40-Mb  (m-get-set L40-Mb  0))  
  (def E40-Mc  (m-get-set L40-Mc  0))  
  (def E40-Pa  (m-get-set L40-Pa  0))
  (def E40-Pc  (m-get-set L40-Pc  0))
  (def E40-S1a (m-get-set L40-S1a 0))
  (def E40-S1b (m-get-set L40-S1b 0)))

(m-get-set 572-Pb.inc  0)
(m-get-set 572-S1.inc  0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET WINDOWED SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#_(defn win-100 [file]
    (p/m-slide-mean file :pi :pi_slide 100))

(defn win-100 [file]
  (p/m-slide-mean file :minfr :minfr_slide 100))

(defn les-win []
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
