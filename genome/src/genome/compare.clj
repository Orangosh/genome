(ns genome.compare
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]
           [genome.database   :as gd ]
           [genome.stats      :as gs ]
           [genome.pop        :as p  ]
           [genome.consvar    :as cv ]
           [genome.dna2aa     :as da ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;
;; For PC
;;(def home "/home/yosh/datafiles/incanted_files/")

;;;;;;;;;;;;;;;;;;;;;;;
;; For Server
(def home "/mnt/data/datafiles/incanted_files/")


(def L05-Pa  (str home "505-Pa.inc" ))
(def L05-M   (str home "505-M.inc"  ))

(def L19-Pb  (str home "519-Pb.inc" ))
(def L19-Pc  (str home "519-Pc.inc" ))
(def L19-Pd  (str home "519-Pd.inc" ))
(def L19-S1a (str home "519-S1a.inc"))

(def L20-Pa  (str home "520-Pa.inc" ))
(def L20-Pb  (str home "520-Pb.inc" ))
(def L20-Pc  (str home "520-Pc.inc" ))
(def L20-S1  (str home "520-S1.inc" )) 
(def L20-S1a (str home "520-S1a.inc"))
  
(def L79-Pa  (str home "579-Pa.inc" ))
(def L79-Pb  (str home "579-Pb.inc" ))
(def L79-M   (str home "579-M.inc"  ))
(def L79-S1a (str home "579-S1a.inc"))
(def L79-S1b (str home "579-S1b.inc"))


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
(def S20-S1a (m-get-set L20-S1a 0)) 
(def S20-S1b (m-get-set L20-S1  0))
  
(def S79-Pa  (m-get-set L79-Pa  0))
(def S79-Pb  (m-get-set L79-Pb  0))
(def S79-M   (m-get-set L79-M   0))
(def S79-S1a (m-get-set L79-S1a 0))
(def S79-S1b (m-get-set L79-S1b 0))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TESTING SNP TREND
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn snp-precent [snp file]
  "snp sould be the string of the keyword"
  (->> file
       (i/add-derived-column
        (keyword (str snp "fq"))
        [(keyword snp) :depth]
        #(if (= %2 0.0)
           0.0
           (/ %1 %2)))))



(defn add-snp-precent [file]
  (->> file
       (snp-precent "T")
       (snp-precent "C")
       (snp-precent "A")
       (snp-precent "G")))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE NEW COLUMN NAME
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn  col-rename [file ser_num]
  "This one is for adding a serial number at the end of a col name"
  (let [dont-change #{:loc  :merlin :gfwd+ :gfwd- :gbwd+ :gbwd-
                      :CDS+ :CDS-   :exon+ :exon-}
        old-cols (apply vector (remove #(contains? dont-change %)
                                       (i/col-names file)))
        new-cols (->> old-cols
                      (map #(keyword (subs (str % ser_num) 1)))
                      (apply vector))]
    (->> new-cols
         (interleave old-cols)
         (apply assoc {}))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;UNITE TWO DATASET AT COMMON SITES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn add-row [file]
  "adds nil as first row because of the join bug"
  (let [ it_be (i/col-names file)]
    (i/conj-rows
     (i/dataset
      it_be
      [(vec (take (count it_be) (repeat nil)))])
     file)))

(defn merge-prep [file]
  "creates a dataset which contains all sites with allele frequency"
  (add-row (add-snp-precent file)))

(defn unite
  ([file1 file2]
   "Adds only file2 rows that have a common :loc value with file1"
   (let [set1   (merge-prep file1)
         set2   (merge-prep file2)
         coled1 (i/rename-cols (col-rename set1 1) set1)
         coled2 (i/rename-cols (col-rename set2 2) set2)]
     (->> coled1
          (i/$join [:loc :loc] coled2)
          (i/$where {:A2 {:$ne nil}}))))
  ([file1 file2 file3]
   (let [set1   (merge-prep file1)
         set2   (merge-prep file2)
         set3   (merge-prep file3)
         coled1 (i/rename-cols (col-rename set1 1) set1)
         coled2 (i/rename-cols (col-rename set2 2) set2)
         coled3 (i/rename-cols (col-rename set3 3) set3)]
     (->> coled1
          (i/$join [:loc :loc] coled2)
          (i/$join [:loc :loc] coled3)
          (i/$where {:A2 {:$ne nil}}))))
  ([file1 file2 file3 file4]
   (let [set1   (merge-prep file1)
         set2   (merge-prep file2)
         set3   (merge-prep file3)
         set4   (merge-prep file4)
         coled1 (i/rename-cols (col-rename set1 1) set1)
         coled2 (i/rename-cols (col-rename set2 2) set2)
         coled3 (i/rename-cols (col-rename set3 3) set3)
         coled4 (i/rename-cols (col-rename set4 4) set4)]
     (->> coled1
          (i/$join [:loc :loc] coled2)
          (i/$join [:loc :loc] coled3)
          (i/$join [:loc :loc] coled4)
          (i/$where {:A2 {:$ne nil}}))))) 
(def p-unite (memoize unite))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; After filtering (i/$where) get genes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#_(defn get-united []
    ;;overtime
    (def a79b79          (p-unite S79-Pa S79-Pb                  ))
    (def a20b20c20       (p-unite S20-Pa S20-Pb  S20-Pc          ))
    (def b19c19d19       (p-unite S19-Pb S19-Pc  S19-Pd          ))
    ;;overtime siblings
    (def Sa79Sb79        (p-unite S79-S1a S79-S1b                ))
    (def Sa20Sb20        (p-unite S20-S1a S20-S1b                ))
    ;;
    (def a5b19a20a79     (p-unite S05-Pa  S19-Pb  S20-Pa S79-Pa  ))
    (def Sa19Sa20M79Sa79 (p-unite S19-S1a S20-S1a S79-M  S79-S1a)))

(defn filtre4 [file] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file
       (i/$where (i/$fn [depth1 minfr1 depth2 minfr2
                         depth3 minfr3 depth4 minfr4]
                        (and (not= nil depth1) (not= nil depth2)
                             (not= nil depth3) (not= nil depth4))))
       (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2
                         majorf+3 minorf+3 majorf+4 minorf+4
                         majorf-1 minorf-1 majorf-2 minorf-2
                         majorf-3 minorf-3 majorf-4 minorf-4]
                        (or (and (not= majorf+1 minorf+1)
                                 (not= majorf+2 minorf+2)
                                 (not= majorf+3 minorf+3)
                                 (not= majorf+4 minorf+4))
                            (and (not= majorf-1 minorf-1)
                                 (not= majorf-2 minorf-2)
                                 (not= majorf-3 minorf-3)
                                 (not= majorf-4 minorf-4)))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3 minfr4]
                        (and (> minfr1 0.003) (> minfr2 0.003)
                             (> minfr3 0.003) (> minfr4 0.003))))
       (i/$where (i/$fn [depth1 depth2 depth3 depth4]
                        (and (> depth1 20.0 ) (> depth2 20.0 )
                             (> depth3 20.0 ) (> depth4 20.0 ))))
       (i/$where (i/$fn [CDS+ CDS-]
                        (or  (not= CDS+ "-" ) (not= CDS- "-"))))
       (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS- :1
             :majorf+1 :majorf+2 :majorf+3 :majorf+4  :2
             :minorf+1 :minorf+2 :minorf+3 :minorf+4  :3
             :majorf-1 :majorf-2 :majorf-3 :majorf-4  :4
             :minorf-1 :minorf-2 :minorf-3 :minorf-4  ])))

(defn filtre3 [file] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file
       (i/$where (i/$fn [depth1 minfr1 depth2 minfr2
                         depth3 minfr3 ]
                        (and (not= nil depth1) (not= nil depth2) (not= nil depth3) )))
       (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2 majorf+3 minorf+3
                         majorf-1 minorf-1 majorf-2 minorf-2 majorf-3 minorf-3]
                        (or (and (not= majorf+1 minorf+1)
                                 (not= majorf+2 minorf+2)
                                 (not= majorf+3 minorf+3))
                            (and (not= majorf-1 minorf-1)
                                 (not= majorf-2 minorf-2)
                                 (not= majorf-3 minorf-3)))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3]
                        (and (> minfr1 0.003) (> minfr2 0.003) (> minfr3 0.003))))
       (i/$where (i/$fn [pi1 pi2 pi3]
                        (and (> pi1 0.0) (> pi2 0.0) (> pi3 0.0))))
       (i/$where (i/$fn [depth1 depth2 depth3]
                        (and (> depth1 20.0 ) (> depth2 20.0 ) (> depth3 20.0 ))))
       (i/$where (i/$fn [CDS+ CDS-] (or  (not= CDS+ "-" ) (not= CDS- "-"))))
       (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS-
             :majorf+1 :majorf+2 :majorf+3
             :minorf+1 :minorf+2 :minorf+3
             :majorf-1 :majorf-2 :majorf-3
             :minorf-1 :minorf-2 :minorf-3])))

(defn filtre2 [file] 
  "A prototype for filther removes nil"
  "gets pos nonsyn, with min allele and depth"
  (->> file
       (i/$where (i/$fn [depth1 minfr1 depth2 minfr2]
                        (and (not= nil depth1) (not= nil depth2))))
       (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2 
                         majorf-1 minorf-1 majorf-2 minorf-2]
                        (or (and (not= majorf+1 minorf+1)
                                 (not= majorf+2 minorf+2))
                            (and (not= majorf-1 minorf-1)
                                 (not= majorf-2 minorf-2)))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3]
                        (and (> minfr1 0.003) (> minfr2 0.003))))
       (i/$where (i/$fn [pi1 pi2]
                        (and (> pi1 0.0) (> pi2 0.0))))
       (i/$where (i/$fn [depth1 depth2 depth3]
                        (and (> depth1 20.0 ) (> depth2 20.0 ))))
       (i/$where (i/$fn [CDS+ CDS-]
                        (or  (not= CDS+ "-" ) (not= CDS- "-"))))
       (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS-
             :majorf+1 :majorf+2 :minorf+1 :minorf+2 
             :majorf-1 :majorf-2 :minorf-1 :minorf-2])))


(defn gene-table[filtre file]
  "returns a table which contains genes and their frequencies given a set and filter"
  (let [g+  (i/$ :gfwd+ (filtre file))
        g-  (i/$ :gfwd+ (filtre file))]
    (i/$order :col-1 :desc
              (#(i/conj-cols (vec (keys %)) (vec (vals %)))
               (frequencies (concat g+ g-))))))

(defn intersect [file1 file2 filtre range-x]
  "returns a set of common genes"
  (let [rngx (min (i/nrow (gene-table filtre file1))
                  (i/nrow (gene-table filtre file2)) range-x)
        tab1 (set (i/$ :col-0 (i/$ (range rngx) :all (gene-table filtre file1))))
        tab2 (set (i/$ :col-0 (i/$ (range rngx) :all (gene-table filtre file2))))]
    (clojure.set/intersection tab1 tab2)))            


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALCULATING MUTATION RATE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn mut-rate [file]
  "input.inc should be a synonimous mutation only"
  (let [com_syn+ (->> file
                      (i/$where {:CDS+ {:$ne "-" }})
                      (i/$where (i/$fn [pi1] (= pi1 0.0 )))
                      (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2] 
                                       (= majorf+1 minorf+1 majorf+2 minorf+2))))
        com_syn- (->> file
                      (i/$where {:CDS- {:$ne "-" }})
                      (i/$where (i/$fn [pi1] (= pi1 0.0 )))
                      (i/$where (i/$fn [majorf-1 minorf-1 majorf-2 minorf-2]
                                       (= majorf-1 minorf-1 majorf-2 minorf-2))))
        mut_syn+ (->> com_syn+
                      (i/$where (i/$fn [maj_p+1 maj_p+2 min_p+2]
                                       (or (not= maj_p+1 maj_p+2)
                                           (not= maj_p+1 min_p+2))))
                      i/nrow)
        mut_syn- (->> com_syn-
                     (i/$where (i/$fn [maj_p-1 maj_p-2 min_p-2]
                                      (or (not= maj_p-1 maj_p-2)
                                          (not= maj_p-1 min_p-2))))
                     i/nrow)
        com_sum  (double (+ (i/nrow com_syn+) (i/nrow com_syn-)))]

    (if (= 0.0 com_sum)
      0.0
      (/ (+ mut_syn+ mut_syn-) com_sum)))) 
    
   
