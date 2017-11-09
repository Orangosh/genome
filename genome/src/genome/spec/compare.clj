(ns genome.compare
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]
           [genome.stats      :as gs ]
           [genome.pop        :as p  ]
           [genome.consvar    :as cv ]
           [genome.dna2aa     :as da ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS- needs one of the specs files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

(defn  c-rename [file ser_num]
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

(defn unite [samples]
  (let [renamed (->> samples
                     (map #(merge-prep %))
                     (map #(i/rename-cols (c-rename %2 %1) %2) (iterate inc 1)))
        enamed  (rest  renamed)]
    (->> renamed
         (reduce #(i/$join [:loc :loc] %1 %2) (map enamed))
         (i/$where {:pi1 {:$ne nil}}))))
    (def p-unite (memoize unite))

       
#_(defn overtime []
    (def a79b79          (p-unite [(samples :S79-Pa)
                                   (samples :S79-Pb)]))
    (def a20b20c20       (p-unite [(samples :S20-Pa)
                                   (samples :S20-Pb)
                                   (samples :S20-Pc)]))
    (def b19c19d19       (p-unite [(samples :S19-Pb)
                                   (samples :S19-Pc)
                                   (samples :S19-Pd)])))

#_(defn get-all-primary []
    (def a5ab79bc20bc19 (p-unite [(samples :S05-Pa)
                                  (samples :S79-Pa)
                                  (samples :S79-Pb)
                                  (samples :S20-Pb)
                                  (samples :S20-Pc)
                                  (samples :S19-Pb)
                                  (samples :S19-Pc)])))

#_(defn overtime-sibs []
    (def Sa79Sb79        (p-unite [S79-S1a
                                   S79-S1b]))
    (def Sa20Sb20        (p-unite [S20-S1a
                                   S20-S1b])))
#_(defn prim-others []
    (def a5b19a20a79     (p-unite [S05-Pa
                                   S19-Pb
                                   S20-Pa
                                   S79-Pa]))
    (def Sa19Sa20M79Sa79 (p-unite [S19-S1a
                                   S20-S1a
                                   S79-M
                                   S79-S1a])))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Post union filtering 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn str>sym [col-str number]
  (apply vector (map #(symbol (str col-str %)) (range 1 (inc number)))))

#_(defn fil-hi-div-nonsyn7 [file & {:keys [filnum dp mf]
                       :or {dp 50
                            mf 0.05}}] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file
       (let [depth   (str>syn "depth"   filnum)
             majorf+ (str>syn "majorf+" filnum)
             majorf- (str>syn "majorf-" filnum)
             minorf+ (str>syn "minorf+" filnum)
             minorf- (str>syn "minorf-" filnum)
             minfr   (str>syn "minfr"   filnum)
             depth   (str>syn "depth"   filnum)]
         (i/$where (i/$fn depth
                          (reduce #(and %1 %2)
                                  (map #(not= nil %)) depth)))
         (i/$where (i/$fn (into [] (concat (majorf+ majorf-
                                            minorf+ minorf-
                                            CDS+ CDS-))) 
                          (or (reduce #(and %1 %2)
                                      (map #(not= %1 %2))
                                      majorf+ minorf+)
                              (reduce #(and %1 %2)
                                      (map #(not= %1 %2))
                                      majorf- minorf-))))
         (i/$where (i/$fn minfr
                          (reduce #(and %1 %2)
                                  (map #(> % mf)) minfr)))
         (i/$where (i/$fn depth
                          (reduce #(and %1 %2)
                                  (map #(> % dp)) minfr)))
         (i/$ (into [] concat
                    [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS-a ]
                    (map #(map #(keyword %) %)
                                   majorf+ majorf- minorf+ minorf- depth))))))))


(defn fil-hi-div-nonsyn7 [file & {:keys [dp mf]
                       :or {dp 50
                            mf 0.05}}] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file 
       (i/$where (i/$fn [depth1 minfr1 depth2 minfr2
                         depth3 minfr3 depth4 minfr4
                         depth5 minfr5 depth6 minfr6
                         depth7 minfr7]
                        (and (not= nil depth1) (not= nil depth2)
                             (not= nil depth3) (not= nil depth4)
                             (not= nil depth5) (not= nil depth6)
                             (not= nil depth7))))
       (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2
                         majorf+3 minorf+3 majorf+4 minorf+4
                         majorf+5 minorf+5 majorf+6 minorf+6
                         majorf+7 minorf+7 CDS+
                         majorf-1 minorf-1 majorf-2 minorf-2
                         majorf-3 minorf-3 majorf-4 minorf-4
                         majorf-5 minorf-5 majorf-6 minorf-6
                         majorf-7 minorf-7 CDS-]
                        (or (and (not= majorf+1 minorf+1)
                                 (not= majorf+2 minorf+2)
                                 (not= majorf+3 minorf+3)
                                 (not= majorf+4 minorf+4)
                                 (not= majorf+5 minorf+5)
                                 (not= majorf+6 minorf+6)
                                 (not= majorf+7 minorf+7) (not= CDS+ "-"))
                            (and (not= majorf-1 minorf-1)
                                 (not= majorf-2 minorf-2)
                                 (not= majorf-3 minorf-3)
                                 (not= majorf-4 minorf-4)
                                 (not= majorf-5 minorf-5)
                                 (not= majorf-6 minorf-6)
                                 (not= majorf-7 minorf-7) (not= CDS- "-")))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3 minfr4 minfr5 minfr6 minfr7]
                        (and (> minfr1 mf) (> minfr2 mf)
                             (> minfr3 mf) (> minfr4 mf)
                             (> minfr5 mf) (> minfr6 mf)
                             (> minfr7 mf))))
       (i/$where (i/$fn [depth1 depth2 depth3 depth4 depth5 depth6 depth7]
                        (and (> depth1 dp ) (> depth2 dp )
                             (> depth3 dp ) (> depth4 dp )
                             (> depth5 dp ) (> depth6 dp )
                             (> depth7 dp ))))
       (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS- 
             :majorf+1 :majorf+2 :majorf+3 :majorf+4 :majorf+5 :majorf+6 :majorf+7
             :minorf+1 :minorf+2 :minorf+3 :minorf+4 :minorf+5 :minorf+6 :minorf+7
             :majorf-1 :majorf-2 :majorf-3 :majorf-4 :majorf-5 :majorf-6 :majorf-7 
             :minorf-1 :minorf-2 :minorf-3 :minorf-4 :mimorf-5 :minorf-6 :minorf-7
             :depth1   :depth2   :depth3   :depth4   :depth5   :depth6   :depth7])))

(defn fil-hi-div-syn7 [file & {:keys [dp mf]
                       :or {dp 50
                            mf 0.05}}] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file 
       (i/$where (i/$fn [depth1 minfr1 depth2 minfr2
                         depth3 minfr3 depth4 minfr4
                         depth5 minfr5 depth6 minfr6
                         depth7 minfr7]
                        (and (not= nil depth1) (not= nil depth2)
                             (not= nil depth3) (not= nil depth4)
                             (not= nil depth5) (not= nil depth6)
                             (not= nil depth7))))
       (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2
                         majorf+3 minorf+3 majorf+4 minorf+4
                         majorf+5 minorf+5 majorf+6 minorf+6
                         majorf+7 minorf+7 CDS+
                         majorf-1 minorf-1 majorf-2 minorf-2
                         majorf-3 minorf-3 majorf-4 minorf-4
                         majorf-5 minorf-5 majorf-6 minorf-6
                         majorf-7 minorf-7 CDS-]
                        (or (and (= majorf+1 minorf+1)
                                 (= majorf+2 minorf+2)
                                 (= majorf+3 minorf+3)
                                 (= majorf+4 minorf+4)
                                 (= majorf+5 minorf+5)
                                 (= majorf+6 minorf+6)
                                 (= majorf+7 minorf+7) (not= CDS+ "-"))
                            (and (= majorf-1 minorf-1)
                                 (= majorf-2 minorf-2)
                                 (= majorf-3 minorf-3)
                                 (= majorf-4 minorf-4)
                                 (= majorf-5 minorf-5)
                                 (= majorf-6 minorf-6)
                                 (= majorf-7 minorf-7) (not= CDS- "-")))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3 minfr4 minfr5 minfr6 minfr7]
                        (and (> minfr1 mf) (> minfr2 mf)
                             (> minfr3 mf) (> minfr4 mf)
                             (> minfr5 mf) (> minfr6 mf)
                             (> minfr7 mf))))
       (i/$where (i/$fn [depth1 depth2 depth3 depth4 depth5 depth6 depth7]
                        (and (> depth1 dp ) (> depth2 dp )
                             (> depth3 dp ) (> depth4 dp )
                             (> depth5 dp ) (> depth6 dp )
                             (> depth7 dp ))))
       (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS- 
             :majorf+1 :majorf+2 :majorf+3 :majorf+4 :majorf+5 :majorf+6 :majorf+7
             :minorf+1 :minorf+2 :minorf+3 :minorf+4 :minorf+5 :minorf+6 :minorf+7
             :majorf-1 :majorf-2 :majorf-3 :majorf-4 :majorf-5 :majorf-6 :majorf-7 
             :minorf-1 :minorf-2 :minorf-3 :minorf-4 :mimorf-5 :minorf-6 :minorf-7
             :depth1   :depth2   :depth3   :depth4   :depth5   :depth6   :depth7])))


(defn fil-low-div7 [file & {:keys [dp mf]
                       :or {dp 50
                            mf 0.0}}] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file 
       (i/$where (i/$fn [depth1 minfr1 depth2 minfr2
                         depth3 minfr3 depth4 minfr4
                         depth5 minfr5 depth6 minfr6
                         depth7 minfr7]
                        (and (not= nil depth1) (not= nil depth2)
                             (not= nil depth3) (not= nil depth4)
                             (not= nil depth5) (not= nil depth6)
                             (not= nil depth7))))
       (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2
                         majorf+3 minorf+3 majorf+4 minorf+4
                         majorf+5 minorf+5 majorf+6 minorf+6
                         majorf+7 minorf+7 CDS+
                         majorf-1 minorf-1 majorf-2 minorf-2
                         majorf-3 minorf-3 majorf-4 minorf-4
                         majorf-5 minorf-5 majorf-6 minorf-6
                         majorf-7 minorf-7 CDS-]
                        (or (and (not= majorf+1 minorf+1)
                                 (not= majorf+2 minorf+2)
                                 (not= majorf+3 minorf+3)
                                 (not= majorf+4 minorf+4)
                                 (not= majorf+5 minorf+5)
                                 (not= majorf+6 minorf+6)
                                 (not= majorf+7 minorf+7) (not= CDS+ "-"))
                            (and (not= majorf-1 minorf-1)
                                 (not= majorf-2 minorf-2)
                                 (not= majorf-3 minorf-3)
                                 (not= majorf-4 minorf-4)
                                 (not= majorf-5 minorf-5)
                                 (not= majorf-6 minorf-6)
                                 (not= majorf-7 minorf-7) (not= CDS- "-")))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3 minfr4 minfr5 minfr6 minfr7]
                        (and (< minfr1 mf) (< minfr2 mf)
                             (< minfr3 mf) (< minfr4 mf)
                             (< minfr5 mf) (< minfr6 mf)
                             (< minfr7 mf))))
       (i/$where (i/$fn [depth1 depth2 depth3 depth4 depth5 depth6 depth7]
                        (and (> depth1 dp ) (> depth2 dp )
                             (> depth3 dp ) (> depth4 dp )
                             (> depth5 dp ) (> depth6 dp )
                             (> depth7 dp ))))
       (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS- 
             :majorf+1 :majorf+2 :majorf+3 :majorf+4 :majorf+5 :majorf+6 :majorf+7
             :minorf+1 :minorf+2 :minorf+3 :minorf+4 :minorf+5 :minorf+6 :minorf+7
             :majorf-1 :majorf-2 :majorf-3 :majorf-4 :majorf-5 :majorf-6 :majorf-7 
             :minorf-1 :minorf-2 :minorf-3 :minorf-4 :mimorf-5 :minorf-6 :minorf-7
             :depth1   :depth2   :depth3   :depth4   :depth5   :depth6   :depth7])))


(defn fil-hi-div-noncoding7 [file & {:keys [dp mf]
                       :or {dp 10
                            mf 0.01}}] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file 
       (i/$where (i/$fn [depth1 minfr1 depth2 minfr2
                         depth3 minfr3 depth4 minfr4
                         depth5 minfr5 depth6 minfr6
                         depth7 minfr7]
                        (and (not= nil depth1) (not= nil depth2)
                             (not= nil depth3) (not= nil depth4)
                             (not= nil depth5) (not= nil depth6)
                             (not= nil depth7))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3 minfr4 minfr5 minfr6 minfr7]
                        (and (> minfr1 mf) (> minfr2 mf)
                             (> minfr3 mf) (> minfr4 mf)
                             (> minfr5 mf) (> minfr6 mf)
                             (> minfr7 mf))))
       (i/$where (i/$fn [depth1 depth2 depth3 depth4 depth5 depth6 depth7]
                        (and (> depth1 dp ) (> depth2 dp )
                             (> depth3 dp ) (> depth4 dp )
                             (> depth5 dp ) (> depth6 dp )
                             (> depth7 dp ))))
       (i/$where (i/$fn [CDS- CDS+]
                        (and (= CDS+ "-") (= CDS- "-"))))
       (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS- 
             :depth1   :depth2   :depth3   :depth4   :depth5   :depth6   :depth7])))



(defn fil-hi-div1 [file & {:keys [dp mf]
                       :or {dp 50
                            mf 0.05}}] 
  "A prototype for filther removes nil, 
gets pos nonsyn, with min allele and depth"
  (->> file 
       (i/$where (i/$fn [depth] (not= nil depth)))
       (i/$where (i/$fn [majorf+ minorf+ CDS+
                         majorf- minorf- CDS-]
                        (or (and (not= majorf+ minorf+)
                                 (not= CDS+ "-"))
                            (and (not= majorf- minorf-)
                                 (not= CDS- "-")))))
       (i/$where (i/$fn [minfr]
                        (and (> minfr mf))))
       (i/$where (i/$fn [depth]
                        (and (> depth dp ))))
       (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS- 
             :majorf+  :majorf- :minorf+ :minorf- :depth])))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;Post filering union and summary
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-genes-table [file filter_name]
  (let [merged (merge-with
                +
                (frequencies (i/$ :gfwd+
                                  (i/$where (i/$fn [CDS+] (not= CDS+ "-"))
                                            (filter_name file))))
                (frequencies (i/$ :gfwd-
                                  (i/$where (i/$fn [CDS-] (not= CDS- "-"))
                                            (filter_name file)))))]
    (->> (i/$order :col-1 :desc (i/conj-cols (keys merged) (vals merged)))
         (i/rename-cols {:col-0 :gene :col-1 :SNPs})
          (i/$where (i/$fn [SNPs] (> SNPs 1)))))) 



























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
                        (and (> minfr1 mf) (> minfr2 mf)
                             (> minfr3 mf) (> minfr4 mf))))
       (i/$where (i/$fn [depth1 depth2 depth3 depth4]
                        (and (> depth1 dp ) (> depth2 dp )
                             (> depth3 dp ) (> depth4 dp ))))
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
                        (and (not= nil depth1)
                             (not= nil depth2)
                             (not= nil depth3) )))
       (i/$where (i/$fn [majorf+1 minorf+1 majorf+2 minorf+2 majorf+3 minorf+3
                         majorf-1 minorf-1 majorf-2 minorf-2 majorf-3 minorf-3]
                        (or (and (not= majorf+1 minorf+1)
                                 (not= majorf+2 minorf+2)
                                 (not= majorf+3 minorf+3))
                            (and (not= majorf-1 minorf-1)
                                 (not= majorf-2 minorf-2)
                                 (not= majorf-3 minorf-3)))))
       (i/$where (i/$fn [minfr1 minfr2 minfr3]
                        (and (> minfr1 mf) (> minfr2 mf) (> minfr3 mf))))
       (i/$where (i/$fn [pi1 pi2 pi3]
                        (and (> pi1 0.0) (> pi2 0.0) (> pi3 0.0))))
       (i/$where (i/$fn [depth1 depth2 depth3]
                        (and (> depth1 dp ) (> depth2 dp ) (> depth3 dp ))))
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
                                 (not= majorf-2 minorf-2))))))
  (i/$where (i/$fn [minfr1 minfr2 minfr3]
                   (and (> minfr1 mf) (> minfr2 mf)))
            (i/$where (i/$fn [pi1 pi2]
                             (and (> pi1 0.0) (> pi2 0.0))))
            (i/$where (i/$fn [CDS+ CDS-]
                             (or  (not= CDS+ "-" ) (not= CDS- "-"))))
            (i/$ [:loc :ref-loc1 :gfwd+ :gfwd- :CDS+ :CDS-
                  :majorf+1 :majorf+2 :minorf+1 :minorf+2 
                  :majorf-1 :majorf-2 :minorf-1 :minorf-2])))






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;MC calculation of Pi
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn circle [[x y]]
  (< (+ (* x x) (* y y)) 1)) 
(defn rand-seq [x] (take x (repeatedly (fn [] (vec (take 2 (repeatedly (fn [] (- (rand 2) 1)))))))))

(defn mcpi [rep]
  (->> (map circle (rand-seq rep))
       frequencies
       ((fn [x] (/ (x true) rep)))
       ((fn [y] (* 4 y))))))))

(defn str>sym-vec [col-str number]
  (->> ()))

