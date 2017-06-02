(ns genome.phylo
  (require [clojure.java.io   :as io ]
           [clojure.string    :as str]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [incanter.charts   :as c  ]
           [genome.analyze     :as da]
           [clojure.data.csv  :as csv]))



(defn get-merlin [col gene file]
  "takes a consensus for a gene in "
  (let [sequence (->> file
                      (i/$where {col {:$eq gene}})
                      (i/$where {:merlin {:$ne "-"}})
                      (i/$ :merlin)
                      (apply str)
                      (partition-all 80)
                      (map #(apply str %)))]
    (spit (str "/mnt/data/datafiles/phylo/input" gene ".fas")
          (str ">merlin\n" (apply str (map #(str % "\n") sequence)))
          :append true)))




(defn get-consensus [col gene sample_name file]
  "takes a consensus for a gene in "
  (let [sequence (->> file
                      (i/$where {col {:$eq gene}})
                      (i/$ :maj+)
                      (apply str)
                      (partition-all 80)
                      (map #(apply str %)))
        ;;       (str ">" sample_name "|gene-" gene "|strand-" )]
        name     (str ">" sample_name)]
    (spit (str "/mnt/data/datafiles/phylo/input/" gene ".fas")
          (str name "\n" (apply str (map #(str % "\n") sequence)))
          :append true)))


(defn get-gene-list [col file]
  (->> file
       (i/$where {col {:$ne "-"}})
       (i/$ col)
       distinct))

(defn map-merlin+ [file]
  (map #(get-merlin :gene+ % file)
       (get-gene-list :gene+ file)))

(defn map-merlin- [file]
  (map #(get-merlin :gene- % file)
       (get-gene-list :gene- file)))


(defn map-gene+ [file]
  (let [sample_name (first (i/$ :r_seq file))]
    (map #(get-consensus :gene+ % sample_name file)
         (get-gene-list :gene+ file))))

(defn map-gene- [file]
  (let [sample_name (first (i/$ :r_seq file))]
    (map #(get-consensus :gene- % sample_name file)
         (get-gene-list :gene- file))))


(defn print-fas+ []
  (map #(map-gene+ %) [S05-M  S05-Pa
                      S19-Pb S19-Pc S19-Pd S19-S1a
                      S20-Pa S20-Pb S20-Pc S20-S1  S20-S1a
                      S79-Pa S79-Pb S79-M  S79-S1a S79-S1b]))
(defn print-fas- []
  (map #(map-gene- %) [S05-M  S05-Pa
                      S19-Pb S19-Pc S19-Pd S19-S1a
                      S20-Pa S20-Pb S20-Pc S20-S1  S20-S1a
                      S79-Pa S79-Pb S79-M  S79-S1a S79-S1b]))

(map-merlin S19-Pb)
(print-fas+)
(print-fas-)
