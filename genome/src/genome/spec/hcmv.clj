(ns genome.spec.hcmv
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [incanter.charts   :as c  ]
           [distributions.core :refer :all]
           [genome.dna2aa     :as da ]
           [genome.stats      :as gs ]           
           [genome.consvar    :as cv ]
           [genome.pop        :as p  ]
           [clojure.data.csv  :as csv]
           [genome.view       :as v  ]
           [clojure.pprint    :as pr ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;
;; For PC
;;(def home "/home/yosh/datafiles/incanted_files/")

;;;;;;;;;;;;;;;;;;;;;;;
 ;; For Server

(defn def-samples [f v & {:keys [depth] :or {depth 0}}]
  (let [home    (fn [x] (str "/mnt/data/hcmv/incanted_files/" x ".inc"))
        samples (zipmap (map keyword v) (map home v))]
    (into {} (for [[k v] samples] [k (f v depth)]))));;

(defn get-set [file cov]
  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))

(defn win-100 [file]
  (p/m-slide-mean file :pi :pi_slide 100))

(def vec-s ["505-Pa"  "505-M"   ;;"519-Pa"
            "519-Pb" "519-Pc" "519-Pd"  "519-S1a"
            "520-Pa" "520-Pb" "520-Pc" "520-S1a" "520-S1b"
            "579-Pa" "579-Pb" "579-M"  "579-S1a" "579-S1b"])
 

#_(def samples (def-samples (fn [x y] (win-100 (m-get-set x y))))) 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn bin [n-bins xs]
  "This function creates bins the output is a map where :map-bin key is an seq of 
to which bin does a value belongs. The :bin-val key describes the size of each bin"
  (let [min-x   (apply min xs)
        max-x   (apply max xs)
        range-x (- max-x min-x)
        bin-fn  (fn [x] (-> x
                            (- min-x)
                            (/ range-x)
                            (* n-bins)
                            (int)
                            (min (dec n-bins))))]
    {:map-bin (map bin-fn xs) :bin-val (/ range-x n-bins)}))

(defn get-inference [file-key shape ci n-bins]
  "Takes one sample and creates a cutoff map where :cutoff is the catoff value 
   assuming a poisson likelyhood and a gamma conjugate prior"
  (let [binned      (->> file-key
                          samples 
                          (i/$ :minfr)
                          (bin n-bins))
        sample-mean (st/mean (:map-bin binned))
        cutoff      (icdf (posterior (:map-bin binned) (poisson :rate) 
                                     (gamma shape (/ shape sample-mean))) ci)]
    {:cutoff cutoff :bin-val (:bin-val binned)}))

(defn cutoff-map [samples & {:keys [shape  ci n-bins]
                             :or   {shape  3
                                    ci     0.95
                                    n-bins 2000}}]
  (i/dataset [:samples :cutoff]
             (seq  (zipmap (keys samples)
                           (map #(let [cut-bin (get-inference % shape ci n-bins)]
                                   (* (cut-bin :cutoff) (cut-bin :bin-val)))
                                (keys samples) )))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;

(defn compare-S19-Pa [samples & sam-keys]
  (let [consensus-gene  (map (fn [file]
                               (->> file
                                    (i/$where (i/$fn [ref-loc]
                                                     (and (< ref-loc 16314)
                                                          (> ref-loc 16250))))
                                    (i/$ :maj+)
                                    (reduce str)))
                             (map #(samples %) sam-keys))
        merlin          (->> (samples (first sam-keys))
                             (i/$where (i/$fn [ref-loc]
                                              (and (< ref-loc 16314)
                                                   (> ref-loc 16250))))
                             (i/$ :ncbi)
                             (reduce str))
        zipped          (zipmap sam-keys consensus-gene)
        conj-mer        (conj {:merlin merlin} zipped)
        home            "/mnt/data/hcmv/S19-Pa.fas"]
    (map #(do (spit home 
                    (str ">" (subs (str (:key %))) "\n")
                    :append true)
              (spit home (:val %)
                    :append true)))))
(defn compare-S19-Pa [samples & sam-keys]
  (let [consensus-gene  (map (fn [file]
                               (->> file
                                    (i/$where (i/$fn [ref-loc]
                                                     (and (< ref-loc 16296)
                                                          (> ref-loc 16250))))
                                    (i/$ :maj+)
                                    (reduce str)))
                             (map #(samples %) sam-keys))
        merlin          (->> (samples (first sam-keys))
                             (i/$where (i/$fn [ref-loc]
                                              (and (< ref-loc 16296)
                                                   (> ref-loc 16250))))
                             (i/$ :ncbi)
                             (reduce str))
        zipped          (zipmap sam-keys consensus-gene)
        conj-mer        (conj {:merlin merlin} zipped)
        home            "/mnt/data/hcmv/S19-Pa.fas"]
    (map #(do (spit home 
                    (str ">" (subs (str (first %)) 1) "\n")
                    :append true)
              (spit home (str (last %) "\n")
                    :append true)) conj-mer)))

(compare-S19-Pa samples :S19-Pa :S19-Pb :S19-Pc :S19-Pd :S19-S1a)                            

                            

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;FOR PCA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn get-inference [file & {:keys [shape rate ci]
                            :or {shape 1 rate 0.1}}]
  (let [sampeled (->> file
                      (i/$ :mf))]
    (icdf (posterion data (poisson :rate) (gamma shape rate))) ci))

(defn le-filter [file & {:keys [m d]
                          :or   {m 0.05
                                 d 35.0}}]
  (->> file
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  depth4
                   depth5  depth6  depth7
                   depth8  depth9
                   depth10 depth11
                   depth12
                   depth13 depth14
                   depth15 depth16]
                  (and (not= nil depth1)
                       (not= nil depth2)  (not= nil depth3)  (not= nil depth4)
                       (not= nil depth5)  (not= nil depth6)  (not= nil depth7)
                       (not= nil depth8)  (not= nil depth9)
                       (not= nil depth10) (not= nil depth11)
                       (not= nil depth12)
                       (not= nil depth13) (not= nil depth14)
                       (not= nil depth15) (not= nil depth16))))
       (i/$where (i/$fn
                  [mf1
                   mf2  mf3  mf4
                   mf5  mf6  mf7
                   mf8  mf9
                   mf10 mf11
                   mf12
                   mf13 mf14
                   mf15 mf16]
                  (or  (> mf1 m)
                       (> mf2 m)  (> mf3 m)  (> mf4 m)
                       (> mf5 m)  (> mf6 m)  (> mf7 m)
                       (> mf8 m)  (> mf9 m)
                       (> mf10 m) (> mf11 m)
                       (> mf12 m)
                       (> mf13 m) (> mf14 m)
                       (> mf15 m) (> mf16 m))))
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  depth4
                   depth5  depth6  depth7
                   depth8  depth9
                   depth10 depth11
                   depth12
                   depth13 depth14
                   depth15 depth16]
                  (and (> depth1  d)
                       (> depth2  d) (> depth3  d) (> depth4  d)
                       (> depth5  d) (> depth6  d) (> depth7  d)
                       (> depth8  d) (> depth9  d)
                       (> depth10 d) (> depth11 d)
                       (> depth12 d)
                       (> depth13 d) (> depth14 d)
                       (> depth15 d) (> depth16 d)
                       )))
       (i/$ [:mf1
             :mf2  :mf3  :mf4
             :mf5  :mf6  :mf7
             :mf8  :mf9
             :mf10 :mf11
             :mf12
             :mf13 :mf14
             :mf15 :mf16
             ])))


(defn CMV-SVD-primary [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 1 2 3 4 5 6 7 8 9] 0 projection)
                        (i/$ [0 1 2 3 4 5 6 7 8 9] 1 projection)
                        :title "CMV primary vs mother sibling"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ (range 9 11)  0 projection)
                      (i/$ (range 9 11)  1 projection))
        (c/add-points (i/$ (range 11 16) 0 projection)
                      (i/$ (range 11 16) 1 projection))
                (i/view))))


(defn CMV-SVD-family [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 9] 0 projection)
                        (i/$ [0 9] 1 projection)
                        :title "CMV families"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [1 2 3 11]  0 projection)
                      (i/$ [1 2 3 11]  1 projection))
        (c/add-points (i/$ [4 5 6 12 13] 0 projection)
                      (i/$ [4 5 6 12 13] 1 projection))
        (c/add-points (i/$ [7 8 10 14 15] 0 projection)
                      (i/$ [7 8 10 14 15] 1 projection))
                (i/view))))

(defn CMV-SVD-depth [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 1 5 8 15] 0 projection)
                        (i/$ [0 1 5 8 15] 1 projection)
                        :title "CMV depth"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [2 4 6 7 14]  0 projection)
                      (i/$ [2 4 6 7 14]  1 projection))
        (c/add-points (i/$ [3 11 12 13 10] 0 projection)
                      (i/$ [3 11 12 13 10] 1 projection))
        (c/add-points (i/$ [9] 0 projection)
                      (i/$ [9] 1 projection))
                (i/view))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Gneral description HCMV
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Gneral descriptive functions

 
(defn cutoff [file]
  (i/nrow (i/$where (i/$fn [minfr pi]
                           (and (> minfr 0.05)
                                (> pi 0.0))) file)))
(defn sum-cov [file]
  (i/sum (i/$ :depth file)))

(defn sum-pi [file]
  (i/sum (i/$ :pi (i/$where (i/$fn [minfr] (> minfr 0.05)) file))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Samples dataset
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn sample-table [samples]
  (->> (i/dataset
        [:sample  :player   :date]
        [[(samples :505-Pa ) "Primary" 1]
         [(samples :505-M  ) "Mother"  1] ;;for graphs remove
         [(samples :519-Pb ) "Primary" 1]
         [(samples :519-Pc ) "Primary" 2]
         [(samples :519-Pd ) "Primary" 3]
         [(samples :519-S1a) "Sibling" 1]
         [(samples :520-Pa ) "Primary" 1]
         [(samples :520-Pb ) "Primary" 2]
         [(samples :520-Pc ) "Primary" 3]
         [(samples :520-S1a) "Sibling" 1]
         [(samples :520-S1b) "Sibling" 2]
         [(samples :579-Pa ) "Primary" 1]
         [(samples :579-Pb ) "Primary" 2]
         [(samples :579-M  ) "Mother"  0]
         [(samples :579-S1a) "Sibling" 1]
         [(samples :579-S1b) "Sibling" 2]])
       (i/add-derived-column
        :name
        [:sample]
        #(first (i/$ :r_seq  %)))
       (i/add-derived-column
        :mean-cov
        [:sample]
        #(format "%.2f" (/ (sum-cov %)
                           (i/nrow  %))))
       (i/add-derived-column
        :cov>20
        [:sample]
        #(i/$where (i/$fn [cov] (> cov 20)) %))
       (i/add-derived-column
        :n-seg
        [:cov>20]
        #(if (> (i/nrow %) 0)
           (format "%.4f" (/ (double (cutoff %)) (i/nrow %))) 0.0))
       (i/add-derived-column
        :nuc-div
        [:cov>20]
        #(if (> (i/nrow %) 0)
           (format "%.4f" (/ (sum-pi %)  (i/nrow %))) 0.0))
       (i/add-derived-column
        :sfs
        [:cov>20]
        #(if (> (i/nrow %) 0) (p/bin-sfs 10 (da/get-synonymous %)) 0))
       (i/$ [:name :player :time-pt :mean-cov :n-seg :nuc-div :sample :cov>20 :sfs])))
(def p-sample-table (memoize sample-table))
;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Simple visualizations
;;; need to fix!!!!!!!!!!!!!!!

(defn show-one [sams [column title result]]
  (i/with-data sams
    (i/view (c/bar-chart :name column
                         :title title
                         :group-by :player
                         :x-label "Sample name"
                         :y-label result
                         :legend true
                         :vertical false))))

(defn show-box [sams [column title result]]
  (-> (c/box-plot column
                  :title title
                  :x-label "Sample name"
                  :y-label result
                  :legend true
                  :series-label "Primary"
                  :data (i/$where {:player {:$eq "Primary"}} sams))
      (c/add-box-plot column
                      :data (i/$where {:player {:$eq "Mother"}} sams)
                      :series-label "Mother")
      (c/add-box-plot column
                      :data (i/$where {:player {:$eq "Sibling"}} sams)
                      :series-label "Sibling")
            (i/view)))


(defn show-common [sams [column title result]]
  (i/with-data  (i/$rollup :mean column :player sams)
    (i/view (c/bar-chart :player column
                         :title title
                         :x-label "Sample name"
                         :y-label result
                         :legend true
                         :vertical false))))

(defn show-all [fnc]
  (let [sams (p-samples)]
    (map #(fnc sams %)

         [[:mean-cov "Mean coverage"                   "mean coverage"        ]
          [:nuc-div  "Nucleotide diversity"            "Nucleotide diversity" ] 
          [:n-seg    "Proportion of segregating sites" "Segregation propotion"]])))
#_(show-all show-one)
#_(show-all show-common)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;ALL SAMPLES descriptive ANALYSIS 

(defn run-all [fnc]
  " a function that accept funcinos returns a valur "
  (let [run_vec (i/$ :cov>20 (p-sample-table))]
    (double (/ (reduce + (map #(fnc %) run_vec))
               (reduce + (map #(i/nrow %) run_vec))))))

(defn stat-all []
  (println "mean coverage for all:           " (run-all sum-cov))
  (println "Nucleotide diversity:            " (run-all sum-pi ))
  (println "Segregating Sites per nucleotid: " (run-all cutoff)))
#_ (stat-all)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;SINGLE SAMPLE descriptive ANALYSIS

(defn summary-sfs [file]
  (let [synonymous (da/get-synonymous file)
        syn-sfs    (p/bin-sfs 10 synonymous)
        mean_cov   (/ (sum-cov file)      (i/nrow file))
        nuc_div    (/ (sum-pi file)       (i/nrow file))
        snp>nucs   (double (/ (cutoff file)       (i/nrow file)))]
    (println "\nSummary statistics for" (first (i/$ :r_seq file))":")
    (println "mean coverage:                    " mean_cov)
    (println "Nucleotide diversity:             " nuc_div)
    (println "Segregating Sites per nucleotid:  " snp>nucs)
    (println "Synonymous site frequency spectra:")
    (println (map first  syn-sfs))
    (println (map second syn-sfs))))

(defn stat-sample [fnc]
  (let [run_vec (i/$ :sample (p-samples))]
    (map #(fnc %) run_vec)))
#_(stat-sample summary-sfs)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Non descriptive single sample
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Functions for understanding minor allele cutoff

(defn fr-dist [dep mfr file] 
  "tool for visualizing minor variant frequencies and FP + Depth"
  (i/$where (i/$fn [depth minfr]
                   (and (> depth dep)
                        (< minfr mfr)))  
            (i/$ [:ref-loc  :gfwd+ :gfwd- 
                  :CDS+     :CDS-  :ref :loc 
                  :depth :T :A  :C :G   :minfr :pi ] 
                 file)))

(defn all-freq-view [file1 file2 dep]
  (-> (c/xy-plot   :loc :minfr
                   :x-label "Position" :y-label "Minor variants frequency"
                   :title (str dep " minimal")
                   :data (fr-dist file1 dep 1.0)) 
      (c/add-lines :loc :minfr
                   :data (fr-dist file2 dep 1.0 ))
      (c/add-lines :loc  :minfr
                   :data (fr-dist file1 dep 0.1))
      (c/add-lines :loc  :minfr
                   :data (fr-dist file2 dep 0.1))
      (c/add-lines :loc  :minfr
                   :data (fr-dist file1 dep 0.003))
      (c/add-lines :loc  :minfr
                   :data (fr-dist file2 dep 0.003))
      (c/add-lines :loc 
                   (map #(/ % 100000) (i/$ :depth file1))
                   :data  file1)
      (c/add-lines :loc 
                   (map #(/ % 10000) (i/$ :depth file2))
                   :data  file2)
      (i/view)))   

(defn poisson-nonfilthered [dep mfr file]
  "get all poisson filtered data which is poistive under certain minor freq (mfr)"
  (i/$where (i/$fn [pi] (> pi 0.0)) (fr-dist file dep mfr)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Different Data slicing methods

(defn pi-chart [file]
  (->> file
       (i/$where (i/$fn [CDS+] (not= CDS+ "-")))
       (i/$ [:ref-loc :gfwd+ :gbwd+ :CDS+ :loc :pi
             :depth :maj+ :min+ :maj_aa+ :min_aa+ 
             :T     :A    :C    :G 
             :orf+ :majorf+ :minorf+])
       (i/$where (i/$fn [pi] (= pi 0.0)))))

(defn nonsyn-chart [file]
  (->> file
       (i/$where (i/$fn [CDS+] (not= CDS+ "-")))
       (i/$ [:ref-loc :gfwd+ :CDS+ :loc 
             :depth :maj+ :min+ :maj_aa+ :min_aa+ 
             :T     :A    :C    :G 
             :orf+ :majorf+ :minorf+])
       (i/$where (i/$fn [majorf+ minorf+] (not= majorf+ minorf+)))))

