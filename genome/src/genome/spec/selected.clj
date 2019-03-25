(ns genome.spec.selected
  (:require [clojure.java.io  :as io]
            [clojure.data.csv :as csv]
            [com.rpl.specter  :as s]
            [clojure.string   :as st]))

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



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Heart of ns

(defn get-filtered-alignments
  "Takes the sequence form first_nucleotide to last_nucleotide 
  sample_name_seq is a vector of  names of the samples you'd like analyse
  map_seq is the actual map of :sample_names and fasta strings"
  [fst_nucleotide lst_nucleotide map_seq sample_name_seq]
  (->> sample_name_seq
       (map #(map_seq %))
       (map #(subs % fst_nucleotide lst_nucleotide))
       (map clojure.string/upper-case)
       (filter #(not= (last %) \-)))) ;Removes last char if it is \- so non left in aligned


(defn aligned?
  "This is the function that is in the heart of the recrussion function (find alinged
  basicaly if all 'priimaries' has a common nucleotide 
  and all 'others' has 'min_true_seqs' with seqs common to 'primaries
  function returns true"
  [lst_nucleotide min_true_samples primaries others map_seq]
  (let [filter-fn     (fn [x] (get-filtered-alignments (dec lst_nucleotide)
                                                      lst_nucleotide
                                                      map_seq
                                                      x))
        concatenated (concat primaries others)
        concat_seqs  (filter-fn concatenated)
        primary_seqs         (let [filtered (filter-fn primaries)]
                               (if (> 2 (count filtered))
                                 false ;;if there are <2 samples that have a value
                                 (apply = filtered))) ;are they equal
        all_seqs (map (fn [y] (let [filtered (filter-fn (concat primaries y))]
                                (if (> 2 (count filtered))
                                  false
                                  (apply not= filtered))))
                      (map vector others))]
    (and primary_seqs
         (<= min_true_samples (count (filter #(= % true) all_seqs)))
         (= (count concatenated) (count concat_seqs)))))  

(defn find-aligned
  [min_true_samples primaries others ref map_seq]
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

(defn get-results
  [min_seq_size min_true_samples all_aligned_seqs]
  (let [filtered (filter #(<= min_seq_size (- (second %) (first %))) all_aligned_seqs)]
    {:results  (vec filtered)
     :stats     {:min_common_nucleotid    min_seq_size
                 :min_non_sharing_samples min_true_samples
                 :number_filtered         (count filtered)}}))


(defn get-all-data
  [ primaries others ref map_seq]
  (vec  (for [z (take 30 (iterate inc 1))
              y (for [x (take (inc (count others)) (iterate inc 0))]
                  {:min_true_samples x :seqs (find-aligned x primaries others ref map_seq)})]
          (get-results z (y :min_true_samples) (y :seqs)))))


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



(def hcmv_data (get-all-data hcmv_primaries hcmv_others hcmv_ref hcmv))
(def hhv6_data (get-all-data hhv6_primaries hhv6_others hhv6_ref hhv6))
(def ebv_data  (get-all-data  ebv_primaries  ebv_others  ebv_ref ebv ))

(def hcmv_csv (map :stats hcmv_data))
(def hhv6_csv (map :stats hhv6_data))
(def  ebv_csv (map :stats ebv_data ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Show distribution of snp's and sort funciton this is for the d3 heatmap in git/visuals/visuals;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
  [path results_vals]
  (map #(print-csv path %) results_vals))


(def hcmv_heatmap1 (get-results-vals 1 hcmv hcmv_ref hcmv_data))
(def hcmv_heatmap2 (get-results-vals 2 hcmv hcmv_ref hcmv_data))
(def hcmv_heatmap3 (get-results-vals 3 hcmv hcmv_ref hcmv_data))
(def  ebv_heatmap1 (get-results-vals 1 ebv   ebv_ref  ebv_data))
(def  ebv_heatmap2 (get-results-vals 2 ebv   ebv_ref  ebv_data))
(def  ebv_heatmap3 (get-results-vals 3 ebv   ebv_ref  ebv_data))
(def hhv6_heatmap1 (get-results-vals 1 hhv6 hhv6_ref hhv6_data))
(def hhv6_heatmap2 (get-results-vals 2 hhv6 hhv6_ref hhv6_data))
(def hhv6_heatmap3 (get-results-vals 3 hhv6 hhv6_ref hhv6_data))


(get-csvs "/mnt/data/primaries_common/circos_data/hcmv/hcmv-heatmap1/hcmv-heatmap1_" hcmv_heatmap1)
(get-csvs "/mnt/data/primaries_common/circos_data/hcmv/hcmv-heatmap2/hcmv-heatmap2_" hcmv_heatmap2)
(get-csvs "/mnt/data/primaries_common/circos_data/hcmv/hcmv-heatmap3/hcmv-heatmap3_" hcmv_heatmap3)
(get-csvs "/mnt/data/primaries_common/circos_data/ebv/ebv-heatmap1/ebv-heatmap1_"    ebv_heatmap1)
(get-csvs "/mnt/data/primaries_common/circos_data/ebv/ebv-heatmap2/ebv-heatmap2_"    ebv_heatmap2)
(get-csvs "/mnt/data/primaries_common/circos_data/ebv/ebv-heatmap3/ebv-heatmap3_"    ebv_heatmap3)
(get-csvs "/mnt/data/primaries_common/circos_data/hhv6//hhv6-heatmap1/hhv6-heatmap1_" hhv6_heatmap1)
(get-csvs "/mnt/data/primaries_common/circos_data/hhv6//hhv6-heatmap2/hhv6-heatmap2_" hhv6_heatmap2)
(get-csvs "/mnt/data/primaries_common/circos_data/hhv6//hhv6-heatmap3/hhv6-heatmap3_" hhv6_heatmap3)


