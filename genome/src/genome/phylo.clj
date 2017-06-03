(ns genome.phylo
  (require [clojure.java.io  :as io ]
           [clojure.string   :as str]
           [incanter.core    :as i  ]
           [incanter.io      :as ii ]
           [incanter.stats   :as st ]
           [incanter.charts  :as c  ]
           [clojure.data.csv :as csv]
           [genome.analyze   :as ana]))

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
;GET TO THE POINT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-merlin [col gene file]
  "takes a consensus for a gene in "
  (let [sequence (->> file
                      (i/$where {col {:$eq gene}})
                      (i/$where {:merlin {:$ne "-"}})
                      (i/$ :merlin)
                      (apply str)
                      (partition-all 80)
                      (map #(apply str %)))]
    (spit (str "/mnt/data/datafiles/phylo/input/" gene ".fas")
          (str ">merlin\n" (apply str (map #(str % "\n") sequence)))
          :append true)))




(defn get-consensus [col gene sample_name file]
  "takes a consensus for a gene in "
  (let [sequence (->> file
                      (i/$where {col {:$eq gene}})
                      (i/$ :maj+)
                      (apply str)
                      (partition-all 80)
                      (map #(apply str %)))
        ;;       (str ">" sample_name "|gene-" gene "|strand-" )]
        name     (str ">" sample_name)]
    (spit (str "/mnt/data/datafiles/phylo/input/" gene ".fas")
          (str name "\n" (apply str (map #(str % "\n") sequence)))
          :append true)))


(defn get-gene-list [col file]
  (->> file
       (i/$where {col {:$ne "-"}})
       (i/$ col)
       distinct))

(defn map-merlin+ [file]
  (map #(get-merlin :gene+ % file)
       (get-gene-list :gene+ file)))

(defn map-merlin- [file]
  (map #(get-merlin :gene- % file)
       (get-gene-list :gene- file)))


(defn map-gene+ [file]
  (let [sample_name (first (i/$ :r_seq file))]
    (map #(get-consensus :gene+ % sample_name file)
         (get-gene-list :gene+ file))))

(defn map-gene- [file]
  (let [sample_name (first (i/$ :r_seq file))]
    (map #(get-consensus :gene- % sample_name file)
         (get-gene-list :gene- file))))
 
(defn print-fas+ []
  (map #(map-gene+ %) [S05-M  S05-Pa
                      S19-Pb S19-Pc S19-Pd S19-S1a
                      S20-Pa S20-Pb S20-Pc S20-S1  S20-S1a
                      S79-Pa S79-Pb S79-M  S79-S1a S79-S1b]))

(defn print-fas- []
  (map #(map-gene- %) [S05-M  S05-Pa
                      S19-Pb S19-Pc S19-Pd S19-S1a
                      S20-Pa S20-Pb S20-Pc S20-S1  S20-S1a
                      S79-Pa S79-Pb S79-M  S79-S1a S79-S1b]))


(map-merlin+ S19-Pb)
(map-merlin- S19-Pb)
(print-fas+)
(print-fas-)

