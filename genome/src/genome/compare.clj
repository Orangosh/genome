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
        [(keyword snp) :cov_p]
        #(if (= %2 0.0)
           0.0
           (/ %1 %2)))))

(defn add-snp-precent [file]
  (->> file
       (snp-precent "Tpois")
       (snp-precent "Cpois")
       (snp-precent "Apois")
       (snp-precent "Gpois")))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE NEW COLUMN NAME
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def sample1 {:Tpois :T1 :Tpois-fq :Tfq1
              :Cpois :C1 :Cpois-fq :Cfq1
              :Apois :A1 :Apois-fq :Afq1
              :Gpois :G1 :Gpois-fq :Gfq1
              :cov_p :cov1 :pi_pois :pi1
              :ref-loc :ref-loc1 :maj_un+ :maj_un+1
              :maj_p+  :maj_p+1  :min_p+  :min_p+1
              :maj_aa+ :maj_aa+1 :min_aa+ :min_aa+1
              :maj_p-  :maj_p-1  :min_p-  :min_p-1
              :maj_aa- :maj_aa-1 :min_aa-  :min_aa-1 
              :majorf+ :majorf+1 :minorf+ :minorf+1
              :majorf- :majorf-1 :minorf- :minorf-1})

(def sample2 {:Tpois :T2 :Tpois-fq :Tfq2
              :Cpois :C2 :Cpois-fq :Cfq2
              :Apois :A2 :Apois-fq :Afq2
              :Gpois :G2 :Gpois-fq :Gfq2
              :cov_p :cov2 :pi_pois :pi2
              :ref-loc :ref-loc2 :maj_un+ :maj_un+2
              :maj_p+  :maj_p+2  :min_p+  :min_p+2
              :maj_aa+ :maj_aa+2 :min_aa+ :min_aa+2
              :maj_p-  :maj_p-2  :min_p-  :min_p-2
              :maj_aa- :maj_aa-2 :min_aa-  :min_aa-2 
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
              :ref-loc1 :cov1     :pi1      :maj_un+1        
              :A1       :Afq1     :T1       :Tfq1
              :G1       :Gfq1     :C1       :Cfq1 
              :maj_p+1  :min_p+1  :maj_aa+1 :min_aa+1 
              :maj_p-1  :min_p-1  :maj_aa-1 :min_aa-1
              :majorf+1 :minorf+1 :majorf-1 :minorf-1
              :ref-loc2 :cov2     :pi2      :maj_un+2 
              :A2       :Afq2     :T2       :Tfq2
              :G2       :Gfq2     :C2       :Cfq2 
              :maj_p+2  :min_p+2  :maj_aa+2 :min_aa+2
              :maj_p-2  :min_p-2  :maj_aa-2 :min_aa-2
              :majorf+2 :minorf+2 :majorf-2 :minorf-2]))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FINAL FUNCTIONS FOR VSRIANTS ALLELE CHANGE AND DIV CHANGE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn nuc-variants [file1 file2]
  "Shows alleles from two samples at same site"
  (->>(create-dataset file1 file2)
      (i/$ [:merlin   :loc  :gene+ :gene-
            :CDS+     :CDS- :exon+ :exon-
            :ref-loc1 :cov1 :pi1  
            :A1   :T1 :G1   :C1 
            :ref-loc2 :cov2 :pi2
            :A2   :T2 :G2   :C2 
            :maj_un+1 :maj_un+2])))

(defn aa-variants [file1 file2]
  "Shows alleles from two samples at same site"
  (->>(create-dataset file1 file2)
      (i/$ [:merlin   :loc      :gene+    :gene-
            :CDS+     :CDS-     :exon+    :exon-
            :ref-loc1 :cov1
            :A1       :T1       :G1       :C1 
            :ref-loc2 :cov2
            :A2       :T2       :G2       :C2 
            :pi1      :pi2      :maj_un+1 :maj_un+2
            :maj_p+1  :maj_p+2  :min_p+1  :min_p+2
            :maj_aa+1 :maj_aa+2 :min_aa+1 :min_aa+2
            :maj_p-1  :maj_p-2  :min_p-1  :min_p-2
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
      ( #(i/conj-rows (i/$ [:cov1 :cov2 :ref-loc :Afq1 :Afq2]  %)
                      (i/$ [:cov1 :cov2 :ref-loc :Tfq1 :Tfq2]  %)
                      (i/$ [:cov1 :cov2 :ref-loc :Cfq1 :Cfq2]  %)
                      (i/$ [:cov1 :cov2 :ref-loc :Gfq1 :Gfq2]  %) ))
      (i/rename-cols {:Afq1 :fq1 :Afq2 :fq2})
      (i/$where (i/$fn [fq1 fq2] (or (> fq1 0.0) (> fq2 0.0)))) 
      (i/$where (i/$fn [fq1 fq2] (or (< fq1 1.0) (< fq2 1.0))))                    
      (i/$order :ref-loc1 :asc)))

(defn diversity-change [file1 file2]
  "compares diversity change between two sites"
  (->>(create-dataset file1 file2)
      (i/$ [:loc :ref-loc1 :cov1 :pi1 :cov2 :pi2])
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
    
       
  
    
