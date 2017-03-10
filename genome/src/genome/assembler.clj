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
        lk (last kmer) ;rk last of kmer
        rbsk (rest (drop-last (r_comp kmer)))
        rfk (first (r_comp kmer))
        rlk (last (r_comp kmer))]
    (if (empty? db1)
      (def dbase (conj db1 {(vec bsk) {fk {lk 1}}}))
      (if (contains? (db1 (vec bsk)) fk) ;if contains kmer -1
        (if (contains? ((db1 (vec bsk)) fk) lk) ; if contalins all kmer
          (def dbase (update-in db1 [(vec bsk) fk lk] inc))   ;Add pos edge
            (def dbase (update-in db1 [(vec bsk) fk] assoc lk 1)))
        (if (contains? (db1 (vec rbsk)) rfk)
          (if (contains? ((db1 (vec rbsk)) rfk) rlk)
              (def dbase (update-in db1 [(vec rbsk) rfk rlk] inc))    ;Add negedge
              (def dbase (update-in db1 [(vec rbsk) rfk] assoc rlk 1)))
          (def dbase (merge-with conj db1 {(vec bsk) {fk {lk 1}}}))))))) ;New node



(defn parsline [rds db2 k]
  ;Recursion of a read needs the exsiting dbase
  (if (> (count rds) k)
    (let [ kmer (take k rds)]
      (create_kmers db2 kmer)
      (recur (rest rds) dbase k))))



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


  
