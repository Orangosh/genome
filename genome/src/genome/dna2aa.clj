(ns genome.dna2aa
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.datasets :as id]
           [incanter.io :as ii ]
           [incanter.charts :as c]
           [incanter.stats :as st]
           [clojure.string :as s]
           [clojure.data.csv :as csv]
           [incanter.distributions :as dst]
           [incanter.zoo :as z]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;COMPLEMENTARY DNA 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;adds column of complementary DNA strand
(defn pos>neg [scanned_column column_name file]
  (let [complementary {"A" "T" "T" "A" "C" "G" "G" "C" "-" "-"}]
    (i/add-column
     column_name
     (->> (i/$ scanned_column file)
          (replace complementary))
     file)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET FORWERD and REVERSE AMINO ACIDS 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;add columns with amino acid after transcription and translation.
(defn nuc>aa [scanned_fwd scanned_rev file]
  (let [DNA>protein  {"TTT" "F" "TTC" "F" "TTA" "L" "TTG" "L"
                      "TCT" "S" "TCC" "S" "TCA" "S" "TCG" "S"
                      "TAT" "Y" "TAC" "Y" "TAA" "$" "TAG" "$"
                      "TGT" "C" "TGC" "C" "TGA" "$" "TGG" "W" 
                      "CTT" "L" "CTC" "L" "CTA" "L" "CTG" "L"
                      "CCT" "P" "CCC" "P" "CCA" "P" "CCG" "P"
                      "CAT" "H" "CAC" "H" "CAA" "Q" "CAG" "Q"
                      "CGT" "R" "CGC" "R" "CGA" "R" "CGG" "R"
                      "ATT" "I" "ATC" "I" "ATA" "I" "ATG" "M"
                      "ACT" "T" "ACC" "T" "ACA" "T" "ACG" "T"
                      "AAT" "N" "AAC" "N" "AAA" "K" "AAG" "K"
                      "AGT" "S" "AGC" "S" "AGA" "R" "AGG" "R"
                      "GTT" "V" "GTC" "V" "GTA" "V" "GTG" "V"
                      "GCT" "A" "GCC" "A" "GCA" "A" "GCG" "A"
                      "GAT" "D" "GAC" "D" "GAA" "E" "GAG" "E"
                      "GGT" "G" "GGC" "G" "GGA" "G" "GGG" "G"
                      "TT-" "-" "TC-" "-" "TA-" "-" "TG-" "-" 
                      "T-T" "-" "T-C" "-" "T-A" "-" "T-G" "-"
                      "-TT" "-" "-TC" "-" "-TA" "-" "-TG" "-"
                      "T--" "-" "-T-" "-" "--T" "-"
                      "CT-" "-" "CC-" "-" "CA-" "-" "CG-" "-" 
                      "C-T" "-" "C-C" "-" "C-A" "-" "C-G" "-"
                      "-CT" "-" "-CC" "-" "-CA" "-" "-CG" "-"
                      "C--" "-" "-C-" "-" "--C" "-"
                      "AT-" "-" "AC-" "-" "AA-" "-" "AG-" "-" 
                      "A-T" "-" "A-C" "-" "A-A" "-" "A-G" "-"
                      "-AT" "-" "-AC" "-" "-AA" "-" "-AG" "-"
                      "A--" "-" "-A-" "-" "--A" "-"
                      "GT-" "-" "GC-" "-" "GA-" "-" "GG-" "-" 
                      "G-T" "-" "G-C" "-" "G-A" "-" "G-G" "-"
                      "-GT" "-" "-GC" "-" "-GA" "-" "-GG" "-"
                      "G--" "-" "-G-" "-" "--G" "-"
                      "---" "-"}
        allele       {:maj_p+   :maj_aa+  :maj_p-   :maj_aa-
                      :min_p+   :min_aa+  :min_p-   :min_aa-}]
    (->>(i/add-column
         (allele scanned_fwd)
         (->> (i/$ scanned_fwd file)
              (partition 3 1)
              (map #(apply str %))
              (conj (take (dec 3) (repeat "-")))
              flatten
              (map DNA>protein)
              (replace {nil "-"}))
         file)
        
        (i/add-column
         (allele scanned_rev)
         (->> (i/$ scanned_rev file)
              (partition 3 1)
              (map #(apply str %))
              (cons (take (dec 3) (repeat "-")))
              flatten
              (map s/reverse)
              (map DNA>protein)
              (replace {nil "-"}))
         )))); <-here file should be

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET READING FRAME
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-orf [CDS]
  "Creats a seq of ORF from CDS column"
  (->> CDS
       (partition-by identity)
       (map #(take (count %) (cycle (range 1 4))))))

(defn add-orfs [file]
  "Adds the ORF seq to the inc file"
  (let [CDS+          (flatten (get-orf          (i/$ :CDS+ file))) 
        CDS- (reverse (flatten (get-orf (reverse (i/$ :CDS- file)))))]
    (->> file
         (i/add-column :orf+ CDS+)
         (i/add-column :orf- CDS-))))

(defn aa-orf [orf amino-seq]
  "Gets the codon sizes from get-orf"
  (map #(repeat %1 %2)
       (->> orf
            (map #(partition-all 3 %))
            (apply concat)
            (map count))
       amino-seq))

(defn add-orf-aa [inc_file]
  "Adds column with ORF for each seq"
  (let [file    (add-orfs inc_file)
        orf+    (get-orf (i/$ :CDS+ file)) 
        orf-    (get-orf (reverse (i/$ :CDS- file)))
        majorf+ (->> file
                     (i/$where { :orf+ {:$eq 1}})
                     (i/$ :maj_aa+))
        minorf+ (->> file
                     (i/$where { :orf+ {:$eq 1}})
                     (i/$ :min_aa+))
        majorf- (->> file
                     (i/$where { :orf- {:$eq 1}})
                     (i/$ :maj_aa-)
                     reverse)
        minorf- (->> file
                     (i/$where { :orf- {:$eq 1}})
                     (i/$ :min_aa-)
                     reverse)]
    (->> file
         (i/add-column :majorf+
                       (flatten (aa-orf orf+ majorf+)))
         (i/add-column :minorf+
                       (flatten (aa-orf orf+ minorf+)))
         (i/add-column :majorf-
                       (flatten (reverse (aa-orf orf- majorf-))))
         (i/add-column :minorf-
                       (flatten (reverse (aa-orf orf- minorf-)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET SYNONYMOUS MUTATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-synonymous [file]
  (i/$order
   [:ref-loc]
   :asc
   (i/conj-rows
    (->> file
         (i/$where {:CDS+ {:$ne "-" }})
         (i/$where (i/$fn [majorf+ minorf+] (= majorf+ minorf+))))
    (->> file
         (i/$where {:CDS- {:$ne "-" }})
         (i/$where (i/$fn [majorf- minorf-] (= majorf- minorf-)))))))


  
