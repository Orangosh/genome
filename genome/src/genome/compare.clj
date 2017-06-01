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
;TESTING SNP TREND
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn snp-precent [snp file]
  "snp sould be the string of the keyword"
  (->> file
       (i/add-derived-column
        (keyword (str snp "-fq"))
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

(def sample1 {:T :T1 :T-fq :Tfq1
              :C :C1 :C-fq :Cfq1
              :A :A1 :A-fq :Afq1
              :G :G1 :G-fq :Gfq1 :minfr   :minfr1
              :depth   :depth1   :pi      :pi1
              :ref-loc :ref-loc1 :maj_un+ :maj_un+1
              :maj+    :maj+1    :min+    :min+1
              :maj_aa+ :maj_aa+1 :min_aa+ :min_aa+1
              :maj-    :maj-1    :min-    :min-1
              :maj_aa- :maj_aa-1 :min_aa- :min_aa-1 
              :majorf+ :majorf+1 :minorf+ :minorf+1
              :majorf- :majorf-1 :minorf- :minorf-1})

(def sample2 {:T :T2 :Tfq :Tfq2
              :C :C2 :C-fq :Cfq2
              :A :A2 :A-fq :Afq2
              :G :G2 :G-fq :Gfq2 :minfr   :minfr2
              :depth   :depth2   :pi      :pi2
              :ref-loc :ref-loc2 :maj_un+ :maj_un+2
              :maj+    :maj+2    :min+    :min+2
              :maj_aa+ :maj_aa+2 :min_aa+ :min_aa+2
              :maj-    :maj-2    :min-    :min-2
              :maj_aa- :maj_aa-2 :min_aa- :min_aa-2 
              :majorf+ :majorf+2 :minorf+ :minorf+2
              :majorf- :majorf-2 :minorf- :minorf-2})


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;UNITE TWO DATASET AT COMMON SITES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn add-row [file]
  (let [ it_be (i/col-names file)]
    (i/conj-rows
     (i/dataset
      it_be
      [(vec (take (count it_be) (repeat nil)))])
     file)))


(defn unite [file1 file2]
  "Adds only file2 rows that have a common :loc value with file1"
  (let [set1   (add-row file1)
        set2   (add-row file2)
        coled1 (i/rename-cols sample1 set1)
        coled2 (i/rename-cols sample2 set2)]
    (i/$where {:A2 {:$ne nil}}
              (i/$join [:loc :loc] coled2 coled1))))
(def p-unite (memoize unite))


(defn create-dataset [file1 file2]
  "creates a dataset which contains all sites with allele frequency"
  (let [snp1   (add-snp-precent file1)
        snp2   (add-snp-precent file2)]
    (->>(p-unite snp1 snp2)
        (i/$ [:merlin   :loc      :gene+    :gene-
              :CDS+     :CDS-     :exon+    :exon-
              :ref-loc1 :depth1   :pi1      :maj_un+1        
              :A1       :Afq1     :T1       :Tfq1
              :G1       :Gfq1     :C1       :Cfq1    :minfr1
              :maj+1    :min+1    :maj_aa+1 :min_aa+1 
              :maj-1    :min-1    :maj_aa-1 :min_aa-1
              :majorf+1 :minorf+1 :majorf-1 :minorf-1
              :ref-loc2 :depth2   :pi2      :maj_un+2 
              :A2       :Afq2     :T2       :Tfq2
              :G2       :Gfq2     :C2       :Cfq2    :minfr2 
              :maj+2    :min+2    :maj_aa+2 :min_aa+2
              :maj-2    :min-2    :maj_aa-2 :min_aa-2
              :majorf+2 :minorf+2 :majorf-2 :minorf-2]))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FINAL FUNCTIONS FOR VSRIANTS ALLELE CHANGE AND DIV CHANGE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn nuc-variants [file1 file2]
  "Shows alleles from two samples at same site"
  (->>(create-dataset file1 file2)
      (i/$ [:merlin   :loc  :gene+ :gene-
            :CDS+     :CDS- :exon+ :exon-
            :ref-loc1 :depth1 :pi1 :minfr1
            :A1   :T1 :G1   :C1 
            :ref-loc2 :depth2 :pi2 :minfr2
            :A2   :T2 :G2   :C2 
            :maj+1    :maj+2])))

(defn aa-variants [file1 file2]
  "Shows alleles from two samples at same site"
  (->>(create-dataset file1 file2)
      (i/$ [:merlin   :loc      :gene+    :gene-
            :CDS+     :CDS-     :exon+    :exon-
            :ref-loc1 :depth1   :minfr1
            :A1       :T1       :G1       :C1 
            :ref-loc2 :depth2   :minfr2
            :A2       :T2       :G2       :C2 
            :pi1      :pi2      :maj_un+1 :maj_un+2
            :maj+1    :maj+2    :min+1    :min+2
            :maj_aa+1 :maj_aa+2 :min_aa+1 :min_aa+2
            :maj-1    :maj-2    :min-1    :min-2
            :maj_aa-1 :maj_aa-2 :min_aa-1 :min_aa-2
            :majorf+1 :majorf+2 :minorf+1 :minorf+2
            :majorf-1 :majorf-2 :minorf-1 :minorf-2])))

(defn allele-change [file1 file2]
  "compares all alleles frequencies from two samples at  same site"
  (->>(create-dataset file1 file2)
      (i/$ [:loc  :gene+    :gene-
            :CDS+     :CDS-     
            :cov1 :Afq1 :Tfq1 :Gfq1 :Cfq1
            :cov2 :Afq2 :Tfq2 :Gfq2 :Cfq2])
      ( #(i/conj-rows (i/$ [:depth1 :depth2 :ref-loc :Afq1 :Afq2]  %)
                      (i/$ [:depth1 :depth2 :ref-loc :Tfq1 :Tfq2]  %)
                      (i/$ [:depth1 :depth2 :ref-loc :Cfq1 :Cfq2]  %)
                      (i/$ [:depth1 :depth2 :ref-loc :Gfq1 :Gfq2]  %) ))
      (i/rename-cols {:Afq1 :fq1 :Afq2 :fq2})
      (i/$where (i/$fn [fq1 fq2] (or (> fq1 0.0) (> fq2 0.0)))) 
      (i/$where (i/$fn [fq1 fq2] (or (< fq1 1.0) (< fq2 1.0))))                    
      (i/$order :ref-loc1 :asc)))

(defn diversity-change [file1 file2]
  "compares diversity change between two sites"
  (->>(create-dataset file1 file2)
      (i/$ [:loc :ref-loc1 :depth1 :pi1 :minfr1 :depth2 :pi2 :minfr2])
      (i/add-derived-column
       :gap
       [:pi2 :pi1]
       #(- %1 %2))))



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
    
       
  
    
