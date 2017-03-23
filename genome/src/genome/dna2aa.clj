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
(defn pos>neg [file scanned_column column_name]
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
(defn nuc>aa [file scanned_fwd scanned_rev]
  (let [DNA>protein {"AAA" "F" "AAG" "F" "AAT" "L" "AAC" "L"
                     "AGA" "S" "AGG" "S" "AGT" "S" "AGC" "S"
                     "ATA" "Y" "ATG" "Y" "ATT" "$" "ATC" "$"
                     "ACA" "C" "ACG" "C" "ACT" "$" "ACC" "W" 
                     "GAA" "L" "GAG" "L" "GAT" "L" "GAC" "L"
                     "GGA" "P" "GGG" "P" "GGT" "P" "GGC" "P"
                     "GTA" "H" "GTG" "H" "GTT" "Q" "GTC" "Q"
                     "GCA" "R" "GCG" "R" "GCT" "R" "GCC" "R"
                     "TAA" "I" "TAG" "I" "TAT" "I" "TAC" "M"
                     "TGA" "T" "TGG" "T" "TGT" "T" "TGC" "T"
                     "TTA" "N" "TTG" "N" "TTT" "K" "TTC" "K"
                     "TCA" "S" "TCG" "S" "TCT" "R" "TCC" "R"
                     "CAA" "V" "CAG" "V" "CAT" "V" "CAC" "V"
                     "CGA" "A" "CGG" "A" "CGT" "A" "CGC" "A"
                     "CTA" "D" "CTG" "D" "CTT" "E" "CTC" "E"
                     "CCA" "G" "CCG" "G" "CCT" "G" "CCC" "G"}]
    (->>(i/add-column
         :aa_fwd
         (->> (i/$ scanned_fwd file)
              (partition 3 1)
              (map #(apply str %))
              (conj (take (dec 3) (repeat "-")))
              flatten
              (map DNA>protein)
              (replace {nil "-"}))
         file)
        
        (i/add-column
         :aa_rev
         (->> (i/$ scanned_rev file)
              (partition 3 1)
              (map #(apply str %))
              (cons (take (dec 3) (repeat "-")))
              flatten
              (map s/reverse)
              (map DNA>protein)
              (replace {nil "-"}))
         )))); <-here file should be
