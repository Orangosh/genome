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
