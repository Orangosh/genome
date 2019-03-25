(ns genome.spec.getseqs
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [incanter.charts   :as c  ]
           [genome.pop        :as p  ]
           [clojure.data.csv  :as csv]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;INCANTER FILES ANALYSIS
;;DEFINE FILE LOCATIONS The "loc" functions creates values 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-locs
  "Creates a map of samples name keywords amd relevant dir/filename.inc"  
  [virus_home virus_name samples_names]
  (let [path "/variants_analysis/incanted_files/"]
    (apply hash-map
     (interleave (map #(keyword %) samples_names)
                 (map #(str virus_home virus_name path % ".inc") samples_names)))))

(def virus_home "/mnt/data/")

(def hcmv_primaries_names
  ["505-Pa"
   "519-Pb"  "519-Pc"  "519-Pd"
   "520-Pa"  "520-Pb"  "520-Pc"
   "579-Pa"  "579-Pb"])

(def hcmv_momsib_names
  ["505-M"
   "519-S1a"
   "520-S1a" "520-S1b"
   "579-M"   "579-S1b" "579-S1a" ])


(def ebv_primaries_names
  ["344-Pb"  
   "525-Pa"  
   "538-Pb"
   "540-Pa"  "540-Pc"])

(def ebv_momsib_names
  ["344-S2a" "344-S2b" "344-S2c" "344-S3"  
   "525-Ma"  "525-Mc"  "525-S1a" "525-S1b" "525-S1c" "525-S1d"
   "525-S2a" "525-S2b" "525-S3a" "525-S3b" "525-S3c"   
   "538-Ma"  "538-S1a"
   "540-Ma"  "540-Mb"  "540-Mc"  "540-S1a" "540-S1b"])

(def hhv6_primaries_names
  ["537-Pa"  "537-Pb"  "537-Pc"  
   "542-Pa"  "542-Pb"  "542-Pc"  
   "543-Pb"  "572-Pa" "572-Pb"])

(def hhv6_momsib_names
  ["537-S2a" "537-S2b" "537-S3"
   "542-S1a" "542-S1b" "542-S1c" "542-Ma"  "542-Mb" 
   "543-S1a" "543-S1b" "572-M" "572-S1" ])

(def hcmv_bm_names
  ["S1"
   "S10" "S11" "S12" "S13" "S14" "S15" "S16" "S17" "S18"
   "S23" "S24" "S25" "S26" "S27" "S28" "S29" "S30"])

(def hcmv_loc
  (get-locs virus_home "hcmv"
            (vec (concat hcmv_primaries_names hcmv_momsib_names))))

(def ebv_loc
  (get-locs virus_home "ebv"
            (vec (concat  ebv_primaries_names  ebv_momsib_names))))

(def hhv6_loc
  (get-locs virus_home "hhv6"
            (vec (concat hhv6_primaries_names hhv6_momsib_names))))
  
(def bm_hcmv_loc
  (get-locs virus_home "bmseq"
            hcmv_bm_names))

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;GET SETS-
;;puts datasets into a "samples" maps wich contains
;;keys of samples name and value as an incnter datasest 


(defn get-set [file cov]
  "creats an incanter dataset with data from adress.csv file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))

(defn get-samples [virus_samples_loc cov]
  (apply hash-map
         (flatten (for [[k v] virus_samples_loc]
                    [k (m-get-set v cov)]))))

(defn get-all-samples
  []
  (def    hcmv_samples (get-samples    hcmv_loc 0))
  (def     ebv_samples (get-samples     ebv_loc 0))
  (def    hhv6_samples (get-samples    hhv6_loc 0))
  #_(def bm_hcmv_samples (get-samples bm_hcmv_loc 0)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;CONSENSUS SEQS ANALYSIS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Load and create data

(defn selected
  "Tests for selected areas in primary infection aligned consensus sequenses "
  [aligned_fasta]
  (->> (with-open [rdr (io/reader aligned_fasta)]
         (reduce conj [] (line-seq rdr)))
       (apply hash-map)))

(def hcmv   (selected "/mnt/data/final_con_seq/HCMV-final.fasta"  ))
(def hhv6   (selected "/mnt/data/final_con_seq/HHV6B-final.fasta" ))
(def ebv    (selected "/mnt/data/final_con_seq/EBV-final.fasta"   ))
(def bmhcmv (selected "/mnt/data/final_con_seq/BMHCMV-final.fasta"))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Define nature of samples

(def hcmv_primaries [">579-Pb" ">579-Pa" ">520-Pb" ">520-Pc"
                     ">519-Pd" ">505-Pa" ">519-Pb" ">520-Pa" ">519-Pc"])
(def hcmv_others    [">520-S1a" ">520-S1b" ">579-S1b" ">519-S1a" ">505-M" ">579-S1a"])
(def hcmv_ref        ">Merlin")
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(def hhv6_primaries [">02-537-Pc"  ">02-542-Pc" ">02-543-Pb" ">02-542-Pa"
                     ">02-572-Pb" ">02-572-Pa" ">02-542-Pb" ">02-537-Pb"  ">02-537-Pa"])
(def hhv6_others    [">02-537-S2b"  ">02-572-S1" ">02-537-S2a"
                     ">02-542-Mb" ">02-542-S1c" ">02-542-S1b" ">02-572-M" ">02-537-S3"
                     ">02-542-S1a"  ">02-542-Ma"  ">02-543-S1b" ">02-543-S1a"])
(def hhv6_ref        ">hhv6b-ref-NC_000898.1")
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(def ebv_primaries [">344-Pb" ">525-Pa" ">540-Pa" ">540-Pc"  ">538-Pb"])
(def ebv_others    [ ">344-S2a" ">344-S2b" ">344-S2c" ">344-S3"  ">525-Ma"  ">525-Mc"
                   ">525-S1a" ">525-S1b" ">525-S1c" ">525-S1d" ">525-S2a" ">525-S2b"
                   ">525-S3a" ">525-S3b" ">525-S3c" ">538-Ma"  ">538-S1a" ">540-Ma"
                   ">540-Mb" ">540-Mc"  ">540-S1a" ">540-S1b" ])
(def ebv_ref       ">ebv-ref-NC_007605.1")

