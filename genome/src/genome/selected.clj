(ns selected.core
  (:require [clojure.java.io  :as io]
            [clojure.data.csv :as csv]
            [com.rpl.specter  :as s]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Load and create data

(defn get-filtered-alignments [fst_nucleotide lst_nucleotide map_seq sample_name_seq]
  (->> sample_name_seq
       (map #(map_seq %))
       (map #(subs % fst_nucleotide lst_nucleotide))
       (map clojure.string/upper-case)
       (filter #(not= (last %) \-))))

(defn aligned? [lst_nucleotide min_true_samples primaries others map_seq]
  (let [primary_seqs         (let [filtered  (get-filtered-alignments (dec lst_nucleotide)
                                                                      lst_nucleotide
                                                                      map_seq
                                                                      primaries)]
                               (if (> 2 (count filtered))
                                 false
                                 (apply = filtered)))
        all_seqs (map (fn [y] (let [filtered (get-filtered-alignments (dec lst_nucleotide)
                                                                      lst_nucleotide
                                                                      map_seq
                                                                      (concat primaries y))]
                                (if (> 2 (count filtered))
                                  false
                                  (apply not= filtered))))
                      (map vector others))]
    (and primary_seqs
         (<= min_true_samples (count (filter #(= % true) all_seqs))))))


(defn find-aligned [min_true_samples primaries others ref map_seq]
  (loop [fst_nucleotide   0
         lst_nucleotide   1
         aligned_seq      []
         all_aligned_seqs []]
    (if (< lst_nucleotide (count (map_seq ref)))
      (if (aligned? lst_nucleotide min_true_samples primaries others map_seq)
        (recur fst_nucleotide
               (inc lst_nucleotide)
               [fst_nucleotide lst_nucleotide]
               all_aligned_seqs)
        (recur lst_nucleotide
               (inc lst_nucleotide)
               []
               (if (empty? aligned_seq)
                 all_aligned_seqs
                 (conj all_aligned_seqs aligned_seq))))
      all_aligned_seqs)))

(defn get-results [min_seq_size min_true_samples all_aligned_seqs]
  (let [filtered (filter #(<= min_seq_size (- (second %) (first %))) all_aligned_seqs)]
    {:results  (vec filtered)
     :stats     {:min_common_nucleotid    min_seq_size
                 :min_non_sharing_samples min_true_samples
                 :number_filtered         (count filtered)}}))


(defn get-all-data [ primaries others ref map_seq]
  (vec  (for [z (take 30 (iterate inc 1))
              y (for [x (take (inc (count others)) (iterate inc 0))]
                  {:min_true_samples x :seqs (find-aligned x primaries others ref map_seq)})]
          (get-results z (y :min_true_samples) (y :seqs)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Sync location with reference locations

(defn sync-seqs [genome_ref]
  (loop [gene_seq         genome_ref
         aligned_loc        1
         aligned_seq      [nil]]
    (if (>= (count gene_seq) 1)
        (if (= (subs gene_seq 0 1) "-")
          (recur (subs gene_seq 1)
                 aligned_loc
                 (conj aligned_seq nil))
          (recur (subs gene_seq 1)
                 (inc  aligned_loc)
                 (conj aligned_seq aligned_loc)))
        aligned_seq)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Given a fst and last nucleotide loc gets a string for each sample
(defn get-str-map [map_seq genome_ref results gap]
  (let [locations  (sync-seqs (map_seq genome_ref))
        string-map (map #(let [[fst lst] results
                               [k v]     %      ]
                           [k (subs v (- fst gap) (+ lst gap))])
                        map_seq)]
    (-> string-map
        flatten
        (#(apply hash-map %))
        (assoc :range [(locations (inc (first results)))
                       (locations (inc (last results)))]))))


(defn final-loc
  [map_seq genome_data genome_ref min_non_sharing_location gap]
  (let [results   ((nth genome_data min_non_sharing_location) :results)]
    (mapv #(get-str-map map_seq genome_ref % gap) results)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;

(defn remove-first-fasta
"remove the frist nucleotide form a map of strings for recrussion"
  [fasta_seq]
  (s/transform [s/MAP-VALS] #(subs % 1) fasta_seq))

(defn get-first-fasta
  "Get the first nucleotide from a map of strings for recurssion"
  [fasta_seq]
  (s/transform [s/MAP-VALS] #(subs % 0 1) fasta_seq))

(defn nuc-to-bul
  [first_fasta]
  (let [first_key (ffirst first_fasta)
        first_val (first_fasta first_key)]
    (zipmap
     (keys first_fasta)
     (map #(= first_val (second %)) first_fasta))))

(defn incorp-first-fasta [snp_distribution first_fasta]
  (let [first_bul   (nuc-to-bul first_fasta)
        neg_bul     (zipmap (keys first_bul) (map false? (vals first_bul)))
        cur_val     (snp_distribution first_bul)
        neg_cur_val (snp_distribution neg_bul)]
    (cond (and (nil? cur_val) (nil? neg_cur_val))      (assoc snp_distribution
                                                              first_bul
                                                              1)
          (and (not= nil cur_val) (= nil neg_cur_val)) (assoc snp_distribution
                                                              first_bul
                                                              (inc cur_val))
          (and (= nil cur_val) (not= nil neg_cur_val)) (assoc snp_distribution
                                                              neg_bul
                                                              (inc neg_cur_val)))))


(defn snp-pattern [fas_seq]
  (loop [fasta_seq         fas_seq
         snp_distribution      {}]
    (if (>= (count fasta_seq) 1)
      (recur (remove-first-fasta fasta_seq)
             (incorp-first-fasta snp_distribution (get-first-fasta fasta_seq)))
      snp_distribution)))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
(defn selected
  "Tests for selected areas in primary infection aligned consensus sequenses "
  [aligned_fasta]
  (->> (with-open [rdr (io/reader aligned_fasta)]
         (reduce conj [] (line-seq rdr)))
       (apply hash-map)))

(defn get-table [virus_specific_data]
  (clojure.pprint/print-table (vec (map #( % :stats)
                                        virus_specific_data))))

(def hcmv_primaries [">579-Pb" ">579-Pa" ">520-Pb" ">520-Pc"
                     ">519-Pd" ">505-Pa" ">519-Pb" ">520-Pa" ">519-Pc"])
(def hcmv_others    [">520-S1a" ">520-S1b" ">579-S1b" ">519-S1a" ">505-M" ">579-S1a"])
(def hcmv_ref        ">Merlin")

(def hhv6_primaries [">02-537-Pc"  ">02-542-Pc" ">02-543-Pb" ">02-542-Pa"
                     ">02-572-Pb" ">02-572-Pa" ">02-542-Pb" ">02-537-Pb"  ">02-537-Pa"])
(def hhv6_others    [">02-537-S2b"  ">02-572-S1" ">02-537-S2a"
                     ">02-542-Mb" ">02-542-S1c" ">02-542-S1b" ">02-572-M" ">02-537-S3"
                     ">02-542-S1a"  ">02-542-Ma"  ">02-543-S1b" ">02-543-S1a"])
(def hhv6_ref        ">hhv6b-ref-NC_000898.1")

(def ebv_primaries [">344-Pb" ">525-Pa" ">540-Pa" ">540-Pc"  ">538-Pb"])
(def ebv_others    [ ">344-S2a" ">344-S2b" ">344-S2c" ">344-S3"  ">525-Ma"  ">525-Mc"
                   ">525-S1a" ">525-S1b" ">525-S1c" ">525-S1d" ">525-S2a" ">525-S2b"
                   ">525-S3a" ">525-S3b" ">525-S3c" ">538-Ma"  ">538-S1a" ">540-Ma"
                   ">540-Mb" ">540-Mc"  ">540-S1a" ">540-S1b" ])
(def ebv_ref       ">ebv-ref-NC_007605.1")

(def hcmv (selected "/mnt/data/hhv6/consensus_analysis/final_consensus_seq/HCMV-final.fasta"))

(def hhv6 (selected "/mnt/data/hhv6/consensus_analysis/final_consensus_seq/HHV6B-final.fasta"))

(def ebv  (selected "/mnt/data/hhv6/consensus_analysis/final_consensus_seq/EBV-final.fasta"))

(def bmhcmv  (selected "/mnt/data/hhv6/consensus_analysis/final_consensus_seq/BMHCMV-final.fasta"))

(def hcmv_data (get-all-data hcmv_primaries hcmv_others hcmv_ref hcmv))
(def hhv6_data (get-all-data hhv6_primaries hhv6_others hhv6_ref hhv6))
(def ebv_data  (get-all-data  ebv_primaries  ebv_others  ebv_ref ebv ))

(def hcmv_csv (map :stats hcmv_data))
(def hhv6_csv (map :stats hhv6_data))
(def  ebv_csv (map :stats ebv_data ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Show distribution of snp's and sort funciton this is for the d3 heatmap in git/visuals/visuals

(def snp_pattern_cmv (snp-pattern hcmv))

(defn sort-snp-pattern [snp_pattern]
  (sort (vals snp_pattern)))

(defn get-identical-P
  "Gets only sereies where primaries has the same nucleotide"
  [snp_pattern primaries]
  (filter (fn [x] (apply = (map (fn [y] ((first x) y)) primaries))) snp_pattern))

(def identical_P (get-identical-P snp_pattern_cmv hcmv_primaries))

(sort-snp-pattern identical_P)

(defn print-aligned-fasta[fst_nucleotide lst_nucleotide map_seq file_name]
  (let [region_seq  (fn [sample_genome]
                      (->> sample_genome
                           (#(subs % fst_nucleotide lst_nucleotide))
                            clojure.string/upper-case))]
    (map #(spit (str file_name ".fas")
                (str (first %) "\n" (region_seq (last %)) "\n")
               :append true)
         map_seq)))

(defn write-csv [path row-data]
  (let [columns [:min_common_nucleotid :min_non_sharing_samples :number_filtered ]
        headers (map name columns)
        rows (mapv #(mapv % columns) row-data)]
    (with-open [file (io/writer path)]
      (csv/write-csv file (cons headers rows)))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Creates csv file for the circos heatmap

(defn get-results-vals
  "gets the"
  [x map_seq genome_ref viral_data]
  (let [loc  (sync-seqs (map_seq genome_ref))
        data (->> viral_data
                  (filter (fn [v] (= x (:min_common_nucleotid (:stats v))))))]
    (for [{:keys [results stats]} data]
      (for [[start end] results]
        {:value (stats :min_non_sharing_samples)
         :start start
         :end end}))))

(defn denest-data
  [results-vals]
  (s/select [s/ALL s/ALL] results-vals))

(defn print-csv
  "Input path and relevent value results from the X-data datasets"
  [path results-vals]
  (let [columns [:value :start :end]
        headers (map name columns)
        rows (mapv #(mapv % columns) results-vals)]
    (with-open [file (io/writer (str path ((first results-vals) :value)))]
      (csv/write-csv file (cons headers rows)))))

(defn get-csvs
  "this one creates a single csv file to each mean common nucleotide"
  [path results-vals]
  (map #(print-csv path %) results-vals))


(def hcmv_heatmap (get-results-vals 1 hcmv hcmv_ref hcmv_data))
(def  ebv_heatmap (get-results-vals 1 ebv   ebv_ref  ebv_data))
(def hhv6_heatmap (get-results-vals 1 hhv6 hhv6_ref hhv6_data))

#_(get-csvs "/mnt/data/hcmv/data/hcmv/hcmv-heatmap" hcmv-heatmap)
#_(get-csvs "/mnt/data/hcmv/data/ebv/ebv-heatmap"    ebv-heatmap)
#_(get-csvs "/mnt/data/hcmv/data/hhv6/hhv6-heatmap" hhv6-heatmap)


