(ns d3.analyze
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
;; For Server

(def home                  "/mnt/data/"  )

 
(defn hcmv-loc []
  (def L05-Pa  (str home "hcmv/variants_analysis/incanted_files/505-Pa.inc" ))
  (def L05-M   (str home "hcmv/variants_analysis/incanted_files/505-M.inc"  ))

  (def L19-Pb  (str home "hcmv/variants_analysis/incanted_files/519-Pb.inc" ))
  (def L19-Pc  (str home "hcmv/variants_analysis/incanted_files/519-Pc.inc" ))
  (def L19-Pd  (str home "hcmv/variants_analysis/incanted_files/519-Pd.inc" ))
  (def L19-S1a (str home "hcmv/variants_analysis/incanted_files/519-S1a.inc"))

  (def L20-Pa  (str home "hcmv/variants_analysis/incanted_files/520-Pa.inc" ))
  (def L20-Pb  (str home "hcmv/variants_analysis/incanted_files/520-Pb.inc" ))
  (def L20-Pc  (str home "hcmv/variants_analysis/incanted_files/520-Pc.inc" ))
  (def L20-S1a (str home "hcmv/variants_analysis/incanted_files/520-S1a.inc" )) 
  (def L20-S1b (str home "hcmv/variants_analysis/incanted_files/520-S1b.inc"))
  
  (def L79-Pa  (str home "hcmv/variants_analysis/incanted_files/579-Pa.inc" ))
  (def L79-Pb  (str home "hcmv/variants_analysis/incanted_files/579-Pb.inc" ))
  (def L79-M   (str home "hcmv/variants_analysis/incanted_files/579-M.inc"  ))
  (def L79-S1a (str home "hcmv/variants_analysis/incanted_files/579-S1a.inc"))
  (def L79-S1b (str home "hcmv/variants_analysis/incanted_files/579-S1b.inc")))



(defn ebv-loc []
  (def L44-Pb  (str home "ebv/variants_analysis/incanted_files/344-Pb.inc" ))   
  (def L44-S2a (str home "ebv/variants_analysis/incanted_files/344-S2a.inc"))  
  (def L44-S2b (str home "ebv/variants_analysis/incanted_files/344-S2b.inc"))  
  (def L44-S2c (str home "ebv/variants_analysis/incanted_files/344-S2c.inc"))   
  (def L44-S3  (str home "ebv/variants_analysis/incanted_files/344-S3.inc" ))
  
  (def L25-Ma  (str home "ebv/variants_analysis/incanted_files/525-Ma.inc" ))  
  (def L25-Mc  (str home "ebv/variants_analysis/incanted_files/525-Mc.inc" ))  
  (def L25-Pa  (str home "ebv/variants_analysis/incanted_files/525-Pa.inc" ))  
  (def L25-S1a (str home "ebv/variants_analysis/incanted_files/525-S1a.inc"))  
  (def L25-S1b (str home "ebv/variants_analysis/incanted_files/525-S1b.inc"))  
  (def L25-S1c (str home "ebv/variants_analysis/incanted_files/525-S1c.inc"))  
  (def L25-S1d (str home "ebv/variants_analysis/incanted_files/525-S1d.inc"))  
  (def L25-S2a (str home "ebv/variants_analysis/incanted_files/525-S2a.inc"))  
  (def L25-S2b (str home "ebv/variants_analysis/incanted_files/525-S2b.inc"))  
  (def L25-S3a (str home "ebv/variants_analysis/incanted_files/525-S3a.inc"))  
  (def L25-S3b (str home "ebv/variants_analysis/incanted_files/525-S3b.inc"))  
  (def L25-S3c (str home "ebv/variants_analysis/incanted_files/525-S3c.inc"))  

  (def L38-Ma  (str home "ebv/variants_analysis/incanted_files/538-Ma.inc" ))   
  (def L38-Pb  (str home "ebv/variants_analysis/incanted_files/538-Pb.inc" ))   
  (def L38-S1a (str home "ebv/variants_analysis/incanted_files/538-S1a.inc"))  

  (def L40-Ma  (str home "ebv/variants_analysis/incanted_files/540-Ma.inc" ))  
  (def L40-Mb  (str home "ebv/variants_analysis/incanted_files/540-Mb.inc" ))  
  (def L40-Mc  (str home "ebv/variants_analysis/incanted_files/540-Mc.inc" ))  
  (def L40-Pa  (str home "ebv/variants_analysis/incanted_files/540-Pa.inc" ))
  (def L40-Pc  (str home "ebv/variants_analysis/incanted_files/540-Pc.inc" ))
  (def L40-S1a (str home "ebv/variants_analysis/incanted_files/540-S1a.inc"))
  (def L40-S1b (str home "ebv/variants_analysis/incanted_files/540-S1b.inc")))

(defn hhv6-loc []
  (def L37-Pa  (str home "hhv6/variants_analysis/incanted_files/537-Pa.inc" ))
  (def L37-S2b (str home "hhv6/variants_analysis/incanted_files/537-S2b.inc"))        
  (def L37-Pb  (str home "hhv6/variants_analysis/incanted_files/537-Pb.inc" ))
  (def L37-S3  (str home "hhv6/variants_analysis/incanted_files/537-S3.inc" ))       
  (def L37-Pc  (str home "hhv6/variants_analysis/incanted_files/537-Pc.inc" ))            
  (def L37-S2a (str home "hhv6/variants_analysis/incanted_files/537-S2a.inc"))

  (def L42-Mb  (str home "hhv6/variants_analysis/incanted_files/542-Mb.inc" ))
  (def L42-Ma  (str home "hhv6/variants_analysis/incanted_files/542-Ma.inc" ))
  (def L42-Pa  (str home "hhv6/variants_analysis/incanted_files/542-Pa.inc" ))
  (def L42-Pb  (str home "hhv6/variants_analysis/incanted_files/542-Pb.inc" ))
  (def L42-Pc  (str home "hhv6/variants_analysis/incanted_files/542-Pc.inc" ))
  (def L42-S1a (str home "hhv6/variants_analysis/incanted_files/542-S1a.inc"))
  (def L42-S1b (str home "hhv6/variants_analysis/incanted_files/542-S1b.inc"))
  (def L42-S1c (str home "hhv6/variants_analysis/incanted_files/542-S1c.inc"))

  (def L43-Pb  (str home "hhv6/variants_analysis/incanted_files/543-Pb.inc" ))
  (def L43-S1a (str home "hhv6/variants_analysis/incanted_files/543-S1a.inc"))  
  (def L43-S1b (str home "hhv6/variants_analysis/incanted_files/543-S1b.inc"))

  (def L72-M   (str home "hhv6/variants_analysis/incanted_files/572-M.inc"  ))
  (def L72-Pa  (str home "hhv6/variants_analysis/incanted_files/572-Pa.inc" )))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn get-set [file cov]
  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))

(defn hcmv-sets []
  (hcmv-loc)
  (def hcmv-samples {:S05-Pa  (m-get-set L05-Pa  0)
                     :S05-M   (m-get-set L05-M   0)
                     
                     :S19-Pb  (m-get-set L19-Pb  0)
                     :S19-Pc  (m-get-set L19-Pc  0)
                     :S19-Pd  (m-get-set L19-Pd  0)
                     :S19-S1a (m-get-set L19-S1a 0)
                     
                     :S20-Pa  (m-get-set L20-Pa  0)
                     :S20-Pb  (m-get-set L20-Pb  0)
                     :S20-Pc  (m-get-set L20-Pc  0)
                     :S20-S1b (m-get-set L20-S1b 0) 
                     :S20-S1a (m-get-set L20-S1a 0)
                     
                     :S79-Pa  (m-get-set L79-Pa  0)
                     :S79-Pb  (m-get-set L79-Pb  0)
                     :S79-M   (m-get-set L79-M   0)
                     :S79-S1a (m-get-set L79-S1a 0)
                     :S79-S1b (m-get-set L79-S1b 0)}))

(defn ebv-sets []
  (ebv-loc)
  (def ebv-samples {:E44-Pb  (m-get-set L44-Pb  0)   
                     :E44-S2a (m-get-set L44-S2a 0)  
                     :E44-S2b (m-get-set L44-S2b 0)  
                     :E44-S2c (m-get-set L44-S2c 0)   
                     :E44-S3  (m-get-set L44-S3  0)
                     
                     :E25-Ma  (m-get-set L25-Ma  0)  
                     :E25-Mc  (m-get-set L25-Mc  0)  
                     :E25-Pa  (m-get-set L25-Pa  0)  
                     :E25-S1a (m-get-set L25-S1a 0)  
                     :E25-S1b (m-get-set L25-S1b 0)  
                     :E25-S1c (m-get-set L25-S1c 0)  
                     :E25-S1d (m-get-set L25-S1d 0)  
                     :E25-S2a (m-get-set L25-S2a 0)  
                     :E25-S2b (m-get-set L25-S2b 0)  
                     :E25-S3a (m-get-set L25-S3a 0)  
                     :E25-S3b (m-get-set L25-S3b 0)  
                     :E25-S3c (m-get-set L25-S3c 0)  
                     
                     :E38-Ma  (m-get-set L38-Ma  0)   
                     :E38-Pb  (m-get-set L38-Pb  0)   
                     :E38-S1a (m-get-set L38-S1a 0)  
                     
                     :E40-Ma  (m-get-set L40-Ma  0)  
                     :E40-Mb  (m-get-set L40-Mb  0)  
                     :E40-Mc  (m-get-set L40-Mc  0)  
                     :E40-Pa  (m-get-set L40-Pa  0)
                     :E40-Pc  (m-get-set L40-Pc  0)
                     :E40-S1a (m-get-set L40-S1a 0)
                     :E40-S1b (m-get-set L40-S1b 0)}))


(defn hhv6-sets []
  (hhv6-loc)
  (def hhv6-samples {:H37-Pa  (m-get-set L37-Pa  0)
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
;PREPARING DATA FOR d3/CIRCOS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
                  :ref-loc :desc
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
        ncol (i/$ :ref-loc file-)]
    (->> file+
         (i/add-column
          :end
          ncol)
         (i/rename-cols {:ref-loc :start
                         :gfwd+ :gfwdp
                         :gfwd- :gfwdn
                         :CDS+  :CDSp
                         :CDS-  :CDSn})
         (i/$ [:r_seq  col   :start :end
               :gfwdp :gfwdn :CDSp  :CDSn]))))
(def m-get-range-set (memoize get-range-set))

(defn circosing-hcmv-ann []
  (let [ An19-Pb (m-get-set L19-Pb 0)] 
    [[(m-get-range-set :gfwd+ true  An19-Pb)
      (str home "hcmv/data/hcmv/genepos.csv")]
     [(m-get-range-set :gfwd- false An19-Pb)
      (str home "hcmv/data/hcmv/geneneg.csv")]
     [(m-get-range-set :CDS+  true  An19-Pb)
      (str home "hcmv/data/hcmv/CDSpos.csv" )]
     [(m-get-range-set :CDS-  false An19-Pb)
      (str home "hcmv/data/hcmv/CDSneg.csv" )]]))

(defn circosing-ebv-ann []
  (let [ An44-Pb (m-get-set L44-Pb 0)] 
    [[(m-get-range-set :gfwd+ true  An44-Pb)
      (str home "hcmv/data/ebv/genepos.csv")]
     [(m-get-range-set :gfwd- false An44-Pb)
      (str home "hcmv/data/ebv/geneneg.csv")]
     [(m-get-range-set :CDS+  true  An44-Pb)
      (str home "hcmv/data/ebv/CDSpos.csv" )]
     [(m-get-range-set :CDS-  false An44-Pb)
      (str home "hcmv/data/ebv/CDSneg.csv" )]]))

(defn circosing-hhv6-ann []
  (let [An37-Pb (m-get-set L37-Pb 0)] 
    [[(m-get-range-set :gfwd+ true  An37-Pb)
      (str home "hcmv/data/hhv6/genepos.csv")]
     [(m-get-range-set :gfwd- false An37-Pb)
      (str home "hcmv/data/hhv6/geneneg.csv")]
     [(m-get-range-set :CDS+  true  An37-Pb)
      (str home "hcmv/data/hhv6/CDSpos.csv" )]
     [(m-get-range-set :CDS-  false An37-Pb)
      (str home "hcmv/data/hhv6/CDSneg.csv" )]]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Diversity

(defn circos-database [file]
  (->> file
       (i/rename-cols {:r_seq   :sample
                       :ref-loc :position
                       :minfr   :value})
       (i/$ [:sample :position :value])))



(defn circosing-hcmv-primary []     
  [[(hcmv-samples :S05-Pa)
    (str home "hcmv/data/hcmv/S05-Pa.csv")]

   [(hcmv-samples :S19-Pb)
    (str home "hcmv/data/hcmv/S19-Pb.csv")]
   [(hcmv-samples :S19-Pc)
    (str home "hcmv/data/hcmv/S19-Pc.csv")]
   [(hcmv-samples :S19-Pd)
    (str home "hcmv/data/hcmv/S19-Pd.csv")]

   [(hcmv-samples :S20-Pa)
    (str home "hcmv/data/hcmv/S20-Pa.csv")]
   [(hcmv-samples :S20-Pb)
    (str home "hcmv/data/hcmv/S20-Pb.csv")]
   [(hcmv-samples :S20-Pc)
    (str home "hcmv/data/hcmv/S20-Pc.csv")]

   [(hcmv-samples :S79-Pa)
    (str home "hcmv/data/hcmv/S79-Pa.csv")]
   [(hcmv-samples :S79-Pb)
    (str home "hcmv/data/hcmv/S79-Pb.csv")]])

(defn circosing-hcmv-momsib []          
  [[(hcmv-samples :S05-M)
    (str home "hcmv/data/hcmv/S05-M.csv"  )]

   [(hcmv-samples :S19-S1a)
    (str home "hcmv/data/hcmv/S19-S1a.csv")]   

   [(hcmv-samples :S20-S1)
    (str home "hcmv/data/hcmv/S20-S1.csv" )]
   [(hcmv-samples :S20-S1a)
    (str home "hcmv/data/hcmv/S20-S1a.csv")]

   [(hcmv-samples :S79-M)
    (str home "hcmv/data/hcmv/S79-M.csv"  )]
   [(hcmv-samples :S79-S1a)
    (str home "hcmv/data/hcmv/S79-S1a.csv")]
   [(hcmv-samples :S79-S1b)
    (str home "hcmv/data/hcmv/S79-S1b.csv")]])



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn circosing-ebv-primary []     
  [[(ebv-samples :E44-Pb)
    (str home "hcmv/data/ebv/E44-Pb.csv")]

   [(ebv-samples :E25-Pa)
    (str home "hcmv/data/ebv/E25-Pa.csv")]

   [(ebv-samples :E38-Pb)
    (str home "hcmv/data/ebv/E38-Pb.csv")]

   [(ebv-samples :E40-Pa)
    (str home "hcmv/data/ebv/E40-Pa.csv")]
   [(ebv-samples :E40-Pc)
    (str home "hcmv/data/ebv/E40-Pc.csv")]])

(defn circosing-ebv-momsib []          
  [[(ebv-samples :E44-S2a)
    (str home "hcmv/data/ebv/E44-S2a.csv")]
   [(ebv-samples :E44-S2b)
    (str home "hcmv/data/ebv/E44-S2b.csv")]
   [(ebv-samples :E44-S2c)
    (str home "hcmv/data/ebv/E44-S2c.csv")]
   [(ebv-samples :E44-S3)
    (str home "hcmv/data/ebv/E44-S3.csv")]
   
   [(ebv-samples :E25-Ma)
    (str home "hcmv/data/ebv/E25-Ma.csv")]
   [(ebv-samples :E25-Mb)
    (str home "hcmv/data/ebv/E25-Mb.csv")]
   [(ebv-samples :E25-S1a)
    (str home "hcmv/data/ebv/E25-S1a.csv")]
   [(ebv-samples :E25-S1b)
    (str home "hcmv/data/ebv/E25-S1b.csv")]
   [(ebv-samples :E25-S1c)
    (str home "hcmv/data/ebv/E25-S1c.csv")]
   [(ebv-samples :E25-S1d)
    (str home "hcmv/data/ebv/E25-S1d.csv")]
   [(ebv-samples :E25-S2a)
    (str home "hcmv/data/ebv/E25-S2a.csv")]
   [(ebv-samples :E25-S2b)
    (str home "hcmv/data/ebv/E25-S2b.csv")]
   [(ebv-samples :E25-S3a)
    (str home "hcmv/data/ebv/E25-S3a.csv")]
   [(ebv-samples :E25-S3b)
    (str home "hcmv/data/ebv/E25-S3b.csv")]
   [(ebv-samples :E25-S3c)
    (str home "hcmv/data/ebv/E25-S3c.csv")]
   
   [(ebv-samples :E38-Ma)
    (str home "hcmv/data/ebv/E38-Ma.csv")]
   [(ebv-samples :E38-S1a)
    (str home "hcmv/data/ebv/E38-S1a.csv")]

   [(ebv-samples :E40-Ma)
    (str home "hcmv/data/ebv/E40-Ma.csv")]
   [(ebv-samples :E40-Mb)
    (str home "hcmv/data/ebv/E40-Mb.csv")]
   [(ebv-samples :E40-Mc)
    (str home "hcmv/data/ebv/E40-Mc.csv")]
   [(ebv-samples :E40-S1a)
    (str home "hcmv/data/ebv/E40-S1a.csv")]
   [(ebv-samples :E40-S1b)
    (str home "hcmv/data/ebv/E40-S1b.csv")]])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn circosing-hhv6-primary []     
  [[(hhv6-samples :H37-Pa)
    (str home "hcmv/data/hhv6/H37-Pa.csv")]
   [(hhv6-samples :H37-Pb)
    (str home "hcmv/data/hhv6/H37-Pb.csv")]
   [(hhv6-samples :H37-Pc)
    (str home "hcmv/data/hhv6/H37-Pc.csv")]

   [(hhv6-samples :H42-Pa)
    (str home "hcmv/data/hhv6/H42-Pa.csv")]
   [(hhv6-samples :H42-Pb)
    (str home "hcmv/data/hhv6/H42-Pb.csv")]
   [(hhv6-samples :H42-Pc)
    (str home "hcmv/data/hhv6/H42-Pc.csv")]

   [(hhv6-samples :H43-Pb)
    (str home "hcmv/data/hhv6/H43-Pb.csv")]

   [(hhv6-samples :H72-Pa)
    (str home "hcmv/data/hhv6/H72-Pa.csv")]
   ])

(defn circosing-hhv6-momsib []          
  [[(hhv6-samples :H37-S2a)
    (str home "hcmv/data/hhv6/H37-S2a.csv")]
   [(hhv6-samples :H37-S2b)
    (str home "hcmv/data/hhv6/H37-S2b.csv")]
   [(hhv6-samples :H37-S3)
    (str home "hcmv/data/hhv6/H37-S3.csv")]

   [(hhv6-samples :H42-Ma)
    (str home "hcmv/data/hhv6/H42-Ma.csv")]
   [(hhv6-samples :H42-Mb)
    (str home "hcmv/data/hhv6/H42-Mb.csv")]
   [(hhv6-samples :H42-S1a)
    (str home "hcmv/data/hhv6/H42-S1a.csv")]
   [(hhv6-samples :H42-S1b)
    (str home "hcmv/data/hhv6/H42-S1b.csv")]
   [(hhv6-samples :H42-S1c)
    (str home "hcmv/data/hhv6/H42-S1c.csv")]

   [(hhv6-samples :H43-S1a)
    (str home "hcmv/data/hhv6/H43-S1a.csv")]
   [(hhv6-samples :H43-S1b)
    (str home "hcmv/data/hhv6/H43-S1b.csv")]
   [(hhv6-samples :H72-M)
    (str home "hcmv/data/hhv6/H72-M.csv")]])



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

(defn map-cir-ann-sets [circosing]
  (map #(cir-ann %) circosing))

#_(map-circos-sets (circosing......))
