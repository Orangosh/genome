(ns genome.spec.compare
  (require [clojure.java.io    :as io ]
           [incanter.core      :as i  ]
           [incanter.datasets  :as id ]
           [incanter.io        :as ii ]
           [incanter.charts    :as c  ]
           [incanter.stats     :as st ]
           [clojure.string     :as s  ]
           [clojure.data.csv   :as csv]
           [genome.spec.getseqs :refer :all ]))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Loading data (see genome/spec/getseqs)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#_(g/get-all-samples)

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
;; UNITE TWO DATASET AT COMMON SITES
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

       
(defn hcmv-overtime []
  (def a79b79           (p-unite [(hcmv_samples :579-Pa)
                                  (hcmv_samples :579-Pb)]))
  (def a20b20c20        (p-unite [(hcmv_samples :520-Pa)
                                  (hcmv_samples :520-Pb)
                                  (hcmv_samples :520-Pc)]))
  (def b19c19d19        (p-unite [(hcmv_samples :519-Pb)
                                  (hcmv_samples :519-Pc)
                                  (hcmv_samples :519-Pd)])))

(defn hcmv-primary []
  (def a5ab79bc20bc19   (p-unite [(hcmv_samples :505-Pa)
                                  (hcmv_samples :579-Pa)
                                  (hcmv_samples :579-Pb)
                                  (hcmv_samples :520-Pb)
                                  (hcmv_samples :520-Pc)
                                  (hcmv_samples :519-Pb)
                                  (hcmv_samples :519-Pc)])))

(defn ebv-overtime []
  (def a40b40           (p-unite [(ebv_samples :540-Pa)
                                  (ebv_samples :540-Pc)])))

(defn ebv-primary []
  (def a25b38b44ac40    (p-unite [(ebv_samples :525-Pa)
                                  (ebv_samples :538-Pb)
                                  (ebv_samples :344-Pb)
                                  (ebv_samples :540-Pa)
                                  (ebv_samples :540-Pc)])))
(defn hhv6-overtime []
  (def a37b37c37        (p-unite [(hhv6_samples :537-Pa)
                                  (hhv6_samples :537-Pb)
                                  (hhv6_samples :537-Pc)]))
  (def a42b42c42        (p-unite [(hhv6_samples :542-Pa)
                                  (hhv6_samples :542-Pb)
                                  (hhv6_samples :542-Pc)])))

(defn hhv6-primary []
  (def abc37abc42b43a72 (p-unite [(hhv6_samples :537-Pa)
                                  (hhv6_samples :537-Pb)
                                  (hhv6_samples :537-Pc)
                                  (hhv6_samples :542-Pa)
                                  (hhv6_samples :542-Pb)
                                  (hhv6_samples :542-Pc)
                                  (hhv6_samples :543-Pb)
                                  (hhv6_samples :572-Pa)])))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Get depth and minfr means
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-hcmv-primary-col
"Calculatse means for deapth and minfr for all hcmv primary infection"
  [file]
  (->> file
       (i/add-derived-column
        :depth_mean
        [:depth1 :depth2 :depth3 :depth4 :depth5 :depth6 :depth7]
        #(/ (+ %1 %2 %3 %4 %5 %6 %7) 7))
       (i/add-derived-column
        :minfr_mean
        [:minfr1 :minfr2 :minfr3 :minfr4 :minfr5 :minfr6 :minfr7]
        #(/ (+ %1 %2 %3 %4 %5 %6 %7) 7))))

(defn get-ebv-primary-col
  "Calculatse means for deapth and minfr for all ebv primary infection"
  [file]
  (->> file
       (i/add-derived-column
        :depth_mean
        [:depth1 :depth2 :depth3 :depth4 :depth5 ]
        #(/ (+ %1 %2 %3 %4 %5) 5))
       (i/add-derived-column
        :minfr_mean
        [:minfr1 :minfr2 :minfr3 :minfr4 :minfr5]
        #(/ (+ %1 %2 %3 %4 %5) 5))))

(defn get-hhv6-primary-col
  "Calculatse means for deapth and minfr for all hhv6 primary infection"
  [file]
  (->> file
       (i/add-derived-column
        :depth_mean
        [:depth1 :depth2 :depth3 :depth4 :depth5 :depth6 :depth7 :depth8]
        #(/ (+ %1 %2 %3 %4 %5 %6 %7 %8) 8))
       (i/add-derived-column
        :minfr_mean
        [:minfr1 :minfr2 :minfr3 :minfr4 :minfr5 :minfr6 :minfr7 :minfr8]
        #(/ (+ %1 %2 %3 %4 %5 %6 %7 %8) 8))))

(defn circos [circosing]
  (let [[file_in file_out] circosing]
    (with-open [f-out (io/writer file_out)]
      (csv/write-csv f-out [(map name (i/col-names file_in))])
      (csv/write-csv f-out (i/to-list file_in)))))


(defn hcmv-mean-primary-depth-minfr
  [min_depth]
  (hcmv-primary)
  (->> a5ab79bc20bc19 
       (get-hcmv-primary-col)
       (i/$ [:loc :ref-loc1 :depth_mean :minfr_mean])
       (i/$where (i/$fn [depth_mean] (> depth_mean min_depth)))))

(defn ebv-mean-primary-depth-minfr
  [min_depth]
  (ebv-primary)
  (->>  a25b38b44ac40  
       (get-ebv-primary-col)
       (i/$ [:loc :ref-loc1 :depth_mean :minfr_mean])
       (i/$where (i/$fn [depth_mean] (> depth_mean min_depth)))))

(defn hhv6-mean-primary-depth-minfr
  [min_depth]
  (hhv6-primary)
  (->> abc37abc42b43a72
       (get-hhv6-primary-col)
       (i/$ [:loc :ref-loc1 :depth_mean :minfr_mean])
       (i/$where (i/$fn [depth_mean] (> depth_mean min_depth)))))

(circos
 [(hcmv-mean-primary-depth-minfr 30)
  "/mnt/data/primaries_common/circos_data/hcmv/hcmv_mean_primary_depth_minfr"])

(circos
 [(ebv-mean-primary-depth-minfr 30)
  "/mnt/data/primaries_common/circos_data/ebv/ebv_mean_primary_depth_minfr"])

(circos
 [(hhv6-mean-primary-depth-minfr 30)
  "/mnt/data/primaries_common/circos_data/hhv6/hhv6_mean_primary_depth_minfr"])

(i/$ (range 40000 60000) [
                          :loc :ref-loc1
                          :depth1 :A1 :T1 :C1 :G1 :maj+1 :min+1
                          :depth2 :A2 :T2 :C2 :G2 :maj+2 :min+2
                          :depth3 :A3 :T3 :C3 :G3 :maj+3 :min+3
                          :depth4 :A4 :T4 :C4 :G4 :maj+4 :min+4
                          :depth5 :A5 :T5 :C5 :G5 :maj+5 :min+5
                          :depth6 :A6 :T6 :C6 :G6 :maj+6 :min+6
                          :depth7 :A7 :T7 :C7 :G7 :maj+7 :min+7
                          ] a5ab79bc20bc19)















;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Old
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;








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

