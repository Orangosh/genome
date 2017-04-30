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
                      "GGT" "G" "GGC" "G" "GGA" "G" "GGG" "G"}
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
