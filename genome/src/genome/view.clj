(ns genome.view
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;VIEW ONE FILE PROPERTIES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn pi-view
  ([file]
  (i/view (c/xy-plot
         :loc
         :pi
         :x-label "Position"
         :y-label "diversity"
         :title "Pi win size"
        ;:legend true
         :data file)))
  ([file cov_div]
  (-> (c/xy-plot
       :loc
       :pi
       :x-label "Position"
       :y-label "diversity"
       :title "02-519Pb AKA''Blip'"
      ;:legend true
       :data file)
      (c/add-lines
       :loc
       (map #(float (/ % cov_div)) (i/$ :depth file)) ;:cov
       :data file)
      (i/view))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;WIEW COMPARRISSONES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn look
  "create a graph of the value comparison"
  ([column file]
   (i/view (c/xy-plot
            :loc
            column
            :x-label "frequency"
            :y-label "Location"
            :title   "Change in nucelotide diversity per site between 2 samples"
;:legend true
            :data file)))

  ([column early>inter early>late]
   (-> (c/xy-plot
        :loc
        column
        :x-label "frequency"
        :y-label "Location"
        :title   "Change in nucelotide diversity per site between 2 samples"
        :data early>inter)
       (c/add-lines
        :loc
        column
        :data early>late)
       (i/view)))

  ([column early>inter early>late early>mom_sib]
    (-> (c/xy-plot
         :loc
         column
         :x-label "frequency"
         :y-label "Location"
         :title   "Change in nucelotide diversity per site between 2 samples"
         :data early>inter)
        (c/add-lines
         :loc
         column
         :data early>late)
        (c/add-lines
         :loc
         column
         :data early>mom_sib)
        (i/view)))

  ([column file1 file2 file3 file4 file5 file6]
    (-> (c/xy-plot
         :loc
         column
         :x-label "frequency"
         :y-label "Location"
         :title   "Change in nucelotide diversity per site between 2 samples"
         :data file1)
        (c/add-lines
         :loc
         column
         :data file2)
        (c/add-lines
         :loc
         column
         :data file3)
        (c/add-lines
         :loc
         column
         :data file4)
        (c/add-lines
         :loc
         column
         :data file5)
        (c/add-lines
         :loc
         column
         :data file6)
        (i/view))))
(defn look
  "create a graph of the value comparison"
  ([column file]
   (i/view (c/xy-plot
            :loc
            column
            :x-label "frequency"
            :y-label "Location"
            :title   "Change in nucelotide diversity per site between 2 samples"
;:legend true
            :data file)))

  ([column early>inter early>late]
   (-> (c/xy-plot
        :loc
        column
        :x-label "frequency"
        :y-label "Location"
        :title   "Change in nucelotide diversity per site between 2 samples"
        :data early>inter)
       (c/add-lines
        :loc
        column
        :data early>late)
       (i/view)))

  ([column early>inter early>late early>mom_sib]
    (-> (c/xy-plot
         :loc
         column
         :x-label "frequency"
         :y-label "Location"
         :title   "Change in nucelotide diversity per site between 2 samples"
         :data early>inter)
        (c/add-lines
         :loc
         column
         :data early>late)
        (c/add-lines
         :loc
         column
         :data early>mom_sib)
        (i/view)))

  ([column file1 file2 file3 file4 file5 file6]
    (-> (c/xy-plot
         :loc
         column
         :x-label "frequency"
         :y-label "Location"
         :title   "Change in nucelotide diversity per site between 2 samples"
         :data file1)
        (c/add-lines
         :loc
         column
         :data file2)
        (c/add-lines
         :loc
         column
         :data file3)
        (c/add-lines
         :loc
         column
         :data file4)
        (c/add-lines
         :loc
         column
         :data file5)
        (c/add-lines
         :loc
         column
         :data file6)
        (i/view))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;save incfiles
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn save-inc [file file-out]
  "/home/yosh/datafiles/incanted_files/SVDs.inc"
    (with-open [f-out (io/writer file-out)]
      (csv/write-csv f-out [(map name (i/col-names file))])
      (csv/write-csv f-out (i/to-list file))))
xo

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;simple trasfer of tables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(save-inc (i/$ [:name :player :time-pt :mean-cov :n-seg :nuc-div] table)  "/mnt/data/hcmv/table-hcmv")

(defn HCMV-scatter [file]
  (let [projection (ii/read-dataset file :header true)]
    (-> (c/scatter-plot (i/$ :mean-cov projection)
                        (i/$ :nuc-div  projection)
                        :title "HCMV preliminary data"
                        :x-label "Coverage"
                        :y-label "Nucleotide diversity")
        (i/view))))

(HCMV-scatter "/home/yosh/data/table-hcmv")
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-mean [min-val max-val file]
  (let [mymean (->> file
                    (i/$where (i/$fn [depth] (and (<= min-val depth)
                                                  (> max-val depth))))
                    (i/$ :minfr)
                    (#(if (empty? %) 0.0 (st/mean %))))]
    [min-val max-val mymean]))

(defn get-bins [bin_size max_depth]
  (let [bin_num (int (/ max_depth bin_size))]
    (partition 2 1 (map #(double (* (/ % bin_num) max_depth))
                        (range (inc bin_num))))))

(defn binned-single-sample [bin_size max_depth file]
  (map #(get-mean (first %) (second %) file)
       (get-bins bin_size max_depth)))

(defn get-all-together [bin_size max_depth file]
  (let [baset (->> (apply i/conj-rows (map vec (get-bins bin_size max_depth)))
                      (i/rename-cols {:col-0 :min-val :col-1 :max-val}))
        mappd (map #(/ % (count samples))
                    (apply map + (map #(map last (binned-single-sample
                                                  bin_size
                                                  max_depth %))
                                      (vals file))))]
    (i/add-column :minfr mappd baset)))
