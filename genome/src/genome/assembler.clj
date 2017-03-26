(ns genome.assembler);a very ugly code just for ideas
(require '[clj-biosequence.core :as bs] ;; for base functionality and fasta
         '[clj-biosequence.blast :as bl] ;; for BLAST functionality
         '[clj-biosequence.fastq :as fq] ;; for fastq functionality
         '[clj-biosequence.tmhmm :as tm] ;; for a wrapper for TMHMM
         '[clj-biosequence.alphabet :as abc]
         '[clojure.java.io :as io]
         '[clojure.edn :as edn])



(defn r_comp [sequence] (map  {\A \T, \T \A, \C \G, \G \C \N \N} (reverse sequence)))



(defn create_kmers [db1 kmer]
  (let [bsk (rest (drop-last kmer)) ;dlk bi-shaved kmer
        fk (first kmer) 


(defn debrujn [ln db3 k lnum readsnum]
  ;parsed over the reads and calls parsline
  (if (= 0 (mod lnum 1000))
    (println lnum " lines out of " (+ lnum (count ln))))
  (if ( or ( >  readsnum lnum) (= readsnum 0))
    (if (empty? ln)
      (println "THE END!")
      (let [x (first ln)]
        (parsline (:sequence x) db3 k)
        (recur (rest ln) dbase k (inc lnum) readsnum)))
    (do
      (println "Graph construction is complete. " (take 5 dbase))
      (def dbg dbase))))



(defn create_contigs [dbgraph nextk onecontig contigs lnum variab]
  (loop [dbgraph dbgraph
         nextk nextk
         onecontig onecontig
         contigs contigs
         lnum lnum]
    (if (= 0 (mod (count dbgraph) 1000))
      (println lnum "nodes;\ndbgraph count: "(count dbgraph)"\n contigs "(count contigs)))
    (let [bskey (rest nextk) fkey (first nextk)]
      (if (empty? dbgraph)
        (do (def ficontigs (conj contigs onecontig))
            (println "Compossition of contigs is complete." (map count ficontigs)))
        (if (empty? (val (first dbgraph)))
          (recur (dissoc dbgraph (key (first dbgraph)))
                 nextk onecontig contigs lnum);)
          (if (contains? (dbgraph bskey) fkey) ;and not contains just on copy         
            (if (= 1 (.count ((dbgraph bskey) fkey)))
              (recur (update-in dbgraph [(vec bskey)] dissoc fkey) 
                     (conj (vec bskey)(key (first ((dbgraph bskey) fkey))))
                     (conj onecontig (key (first ((dbgraph bskey) fkey))))
                     contigs (inc lnum));)
              (if (= variab "true")
                (recur (update-in dbgraph [(vec bskey)] dissoc fkey) 
                       (conj (vec bskey) (first (keep #(when (= (val %)
                             (apply max (vals ((dbgraph bskey) fkey))))
                             (key %)) ((dbgraph bskey) fkey))))
                       (conj onecontig (first (keep #(when (= (val %)
                             (apply max (vals ((dbgraph bskey) fkey))))
                             (key %)) ((dbgraph bskey) fkey))))
                       contigs (inc lnum))
                (recur (update-in dbgraph [(vec bskey)] dissoc fkey)
                       (vec (cons (key (first (val (first dbgraph))))
                                  (key (first dbgraph))))
                       (vec (cons (key (first (val (first dbgraph))))
                                  (key (first dbgraph))))
                       (conj contigs onecontig) (inc lnum))))
            (recur dbgraph 
                   (vec (cons (key (first (val (first dbgraph))))
                              (key (first dbgraph))))
                   (vec (cons (key (first (val (first dbgraph))))
                              (key (first dbgraph))))
                   (conj contigs onecontig) lnum)))))))



(defn merge_contigs [mcon k]
  (loop [mconx mcon
         mcony mcon
         mconz mcon
         k k
         ln 1]
    (if (= 1 (.count  mconx))
      (do (def merged_contigs mcony)
          (println (count mconx) "\n mconx" mconx "\n mconz "(map count mconz))
          (spit "/home/yosh/Software/clojure/merged_contigs.dat"
                (prn-str merged_contigs)))
      (let [x (first mconx) y (first mcony)]
        (if (< 1 (.count mcony))
          (if ( and (= (take-last (- k 1) x)  (take (- k 1) y))
               (not (= (take-last (- k 1) x)  (take (- k 1) x))))
            (do (defn adjust [xyz]
                (vec (filter #(not (= y %))
                     (replace {x (vec (flatten (conj x (drop (- k 1) y))))} xyz))))
              (recur (adjust mconx) (adjust mcony) (adjust mconz) k (inc ln)))
            (recur mconx (rest mcony) mconz k (inc ln)));)
          (recur (rest mconx) mconz mconz k (inc ln)))))));)

(defn -main [file strk readsnumstring variability]
  ;defines k (preferably odd number) and the fastq file we want to evaluate
  (println "Please wait while we create a de bruijn graph....")
  (with-open [r (bs/bs-reader (fq/init-fastq-file file))]
    (def k (read-string strk))
    (def readsnum (read-string readsnumstring))
    (debrujn (bs/biosequence-seq r) {} k 1 readsnum))
  (println "Creating contigs")
  (if (= variability "true")
    (println "Getting contigs that include the most common kmers")
    (println "Getting contigs that include only one kmer at each possition"))
  (trampoline create_contigs dbg
              (vec (cons (key (first (val (first dbg)))) (key (first dbg))))
              (vec (cons (key (first (val (first dbg)))) (key (first dbg))))
              [] 1 variability) 
  (trampoline merge_contigs ficontigs k))


  
;transition maps for RNA and rev DNA to Protein
(def RNA>protein {"UUU" "F" "UUC" "F" "UUA" "L" "UUG" "L"
                  "UCU" "S" "UCC" "S" "UCA" "S" "UCG" "S"
                  "UAU" "Y" "UAC" "Y" "UAA" "$" "UAG" "$"
                  "UGU" "C" "UGC" "C" "UGA" "$" "UGG" "W" 
                  "CUU" "L" "CUC" "L" "CUA" "L" "CUG" "L"
                  "CCU" "P" "CCC" "P" "CCA" "P" "CCG" "P"
                  "CAU" "H" "CAC" "H" "CAA" "Q" "CAG" "Q"
                  "CGU" "R" "CGC" "R" "CGA" "R" "CGG" "R"
                  "AUU" "I" "AUC" "I" "AUA" "I" "AUG" "M"
                  "ACU" "T" "ACC" "T" "ACA" "T" "ACG" "T"
                  "AAU" "N" "AAC" "N" "AAA" "K" "AAG" "K"
                  "AGU" "S" "AGC" "S" "AGA" "R" "AGG" "R"
                  "GUU" "V" "GUC" "V" "GUA" "V" "GUG" "V"
                  "GCU" "A" "GCC" "A" "GCA" "A" "GCG" "A"
                  "GAU" "D" "GAC" "D" "GAA" "E" "GAG" "E"
                  "GGU" "G" "GGC" "G" "GGA" "G" "GGG" "G"})

(def revDNA>protein {"TTT" "F" "TTC" "F" "TTA" "L" "TTG" "L"
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


(ns genome.annotate
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.stats :as st]
           [clojure.xml :as xml]
           [clojure.zip :as zip]
           [clojure.data.xml :as cx]
           [clojure.data.zip.xml :as dzx]))

;; wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_006273.2&retmode=xml"
;; mv efetch.fcgi?db=nuccore\&id=NC_006273.2\&retmode=xml merlin.xm

(defn get-ids [zipper]
  "Extract specific elements from an XML document"
  (dzx/xml-> zipper
             :IdList
             :Id
             dzx/text))

(defn get-affiliations [zipper]
  "Extract affiliations from PubMed abstracts"
  (map (fn [a b c d] (vector a b c d))
       (dzx/xml-> zipper
                  :GBSeq
                  :GBSeq_feature-table
                  :GBFeature
                  :GBFeature_intervals
                  :GBInterval
                  :GBInterval_from                  
                  dzx/text)
       (dzx/xml-> zipper
                  :GBSeq
                  :GBSeq_feature-table
                  :GBFeature
                  :GBFeature_intervals
                  :GBInterval
                  :GBInterval_to                  
                  dzx/text)
       (dzx/xml-> zipper
                  :GBSeq
                  :GBSeq_feature-table
                  :GBFeature
                  :GBFeature_quals
                  :GBQualifier
                  :GBQualifier_name
                  dzx/text)
       (dzx/xml-> zipper
                  :GBSeq
                  :GBSeq_feature-table
                  :GBFeature
                  :GBFeature_quals
                  :GBQualifier
                  :GBQualifier_value
                  dzx/text)))


(def data "/home/yosh/datafiles/merlin.xml")

(println (get-ids
      (zip/xml-zip
       (xml/parse data))))

(def dat1 (get-affiliations (->> data
                                 xml/parse
                                 zip/xml-zip
                                 zip/down))



(def col_names {:GBFeature_key :feature
                :GBInterval_from :starts
                :GBInterval_to :ends})

(def qualifiers {:GBQualifier_name :description
                 :GBQualifier_value :specification})

(def replace_col_names {:col-0 :feature
                        :col-1 :starts
                        :col-2 :ends})

(def replace_qualifiers {:col-0 :description
                         :col-1 :specification})

(def data "/home/yosh/datafiles/merlin.xml")

(defn str>double [s] ;; ;returnes the first string number as an integer
  (Double. (re-find  #"\d+" s )))

(defn parseXML [ file col_name]
  (vec (for [x (xml-seq (xml/parse( java.io.File. file)))
             :when (= col_name (:tag x))]
         (first (:content x)))))

(defn annotate [data col_names col_replace]
  (->> (map #(parseXML data %) (vec (keys col_names)))
       (apply i/conj-cols)
       (i/rename-cols col_replace)))

(def annotation (annotate data col_names replace_col_names))
(def qual_inc (annotate data qualifiers replace_qualifiers))

(distinct (i/$ :description qual_inc))

(def ld (i/$where {:feature {:$eq "gene" }} annotation))
(def ld (i/$where {:description {:$eq "gene" }} qual_inc))

(defn load-xml-data [xml-file first-data next-data]
  (let [data-map (fn [node]
                   [(:tag node) (first (:content node))])]
    (->>
     (xml/parse xml-file)
     zip/xml-zip
     first-data
     (iterate next-data)
     (take-while #(not (nil? %)))
     (map zip/children)
     (map #(mapcat data-map %))
     (map #(apply array-map %))
     i/to-dataset)))




