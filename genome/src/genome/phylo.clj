(ns genome.phylo
  (require [clojure.java.io  :as io ]
           [clojure.string   :as str]
           [incanter.core    :as i  ]
           [incanter.io      :as ii ]
           [incanter.stats   :as st ]
           [incanter.charts  :as c  ]
           [clojure.data.csv :as csv]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;
;; For PC

(def home                  "/home/yosh/datafiles/")
(def input_file  (str home "genes/merlin.gff3"   ))
(def output_file (str home "genes/merlin.inc"    ))

;;;;;;;;;;;;;;;;;;;;;;;
;; For Server

;;(def home                  "/mnt/data/datafiles/"  )
;;(def input_file  (str home "concensus/merlin.gff3"))
;;(def output_file (str home "concensus/refset.inc" ))


 
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
;GET A TABLE WITH GENE START END FROM gff3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn pairs [type col]
  "parses through the ID line of gff3, creates a map returns a value"
  (->> col
       (re-seq  #"[^;]*")
       (map #(re-seq #"[^=]*" %))
       (map (fn [x] (filter #(not= "" %) x)))
       (filter #(= 2 (count %)))
       flatten
       (apply hash-map)
       (#(% type))))

(defn gene-list [file]
  "/home/yosh/datafiles/genes/merlin.gff3"
  (->> (ii/read-dataset file :delim \tab)
       (i/rename-cols {:col2 :type :col3 :pt1 :col4 :pt2
                       :col6 :strand   :col8 :attributes})
       (i/$where {:type {:eq "gene"}})
       (i/add-derived-column
        :gene
        [:attributes]
        #(pairs "gene" %))(i/$  [:gene :pt1 :pt2])
       (i/to-vect)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;PRINT THE GENE INTO FASTA FILE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-gene [[gene pt1 pt2] sample_name col file]
  (let [sequence (->> file
                      (i/$where (i/$fn [ref-loc]
                                       (and (<= pt1 ref-loc)
                                            (>= pt2 ref-loc))))
                      (i/$ col)
                      (remove #(= "-" %)) 
                      (apply str)
                      (partition-all 80)
                      (map #(apply str %)) )
        name     (str ">" sample_name)]
    (spit (str "/home/yosh/datafiles/trees/treesfromcon/input/" gene ".fas")
          (str name "\n" (apply str (map #(str % "\n") sequence)))
          :append true)))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn fasta>inc [file]
  "/home/yosh/datafiles/trees/HCMVfortree.fasta"
  (let [fasta (with-open [rdr (io/reader
                   file)]
                (->>(line-seq rdr)
                    doall
                    (partition 2)))
        seq-col  (->> fasta
                   (map second)
                   (map #(vec (str/split % #"")))
                   (apply i/conj-cols))
        name-col (->> fasta
                       (map #(subs (first %) 1))
                       vec)
        locless  (i/rename-cols
                  (apply assoc {} (interleave (i/col-names seq-col) name-col))
                  seq-col)]
    (i/add-column :ref-loc (vec (range 1 (inc (i/nrow locless)))) locless)))



(defn inc>gene [gff-list file col]
  (let [sample_name  (subs (str col) 1)]
    (map #(get-gene % sample_name col file)
         gff-list)))


(defn gene>print [gff file]
  "gff3 /home/yosh/datafiles/genes/merlin.gff3"
  "file /home/yosh/datafiles/trees/HCMVfortree.fasta"
  (let [gff-list (gene-list gff)
        inc-file (fasta>inc file)]
    (map #(inc>gene gff-list inc-file %)
         (i/col-names inc-file))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn map-gene [gff [col file]]
  (let [ sample_name (if (= col :merlin) "merlin" (first (i/$ :r_seq file)))]
    (map #(get-gene % sample_name col file)
         (gene-list gff))))

(defn print-fas []
  (map #(map-gene "/home/yosh/datafiles/genes/merlin.gff3" %)
       [[:merlin S19-Pb]
        [:maj+ S05-M ]  [:maj+ S05-Pa ] [:maj+ S19-Pb ] [:maj+ S19-Pc ]
        [:maj+ S19-Pd]  [:maj+ S19-S1a] [:maj+ S20-Pa ] [:maj+ S20-Pb ]
        [:maj+ S20-Pc]  [:maj+ S20-S1 ] [:maj+ S20-S1a] [:maj+ S79-Pa ]
        [:maj+ S79-Pb]  [:maj+ S79-M  ] [:maj+ S79-S1a] [:maj+ S79-S1b]]))
