(ns genome.annotate
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.stats :as st]
           [incanter.io :as ii ]
           [clojure.xml :as xml]))

;; wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_006273.2&retmode=xml"
;; mv efetch.fcgi?db=nuccore\&id=NC_006273.2\&retmode=xml merlin.xm

(def data "/home/yosh/datafiles/merlin.xml")

(def col_names [:GBFeature_key :GBInterval_from :GBInterval_to])

(def replace_col_names {:col-0 :feature
                        :col-1 :starts
                        :col-2 :ends})

(defn parseXML [ file col_name]
  ""
  (vec (for [x (xml-seq (xml/parse( java.io.File. file)))
             :when (= col_name (:tag x))]
         (first (:content x)))))

(defn annotate [data col_names col_replace]
  ""
  (->> (map #(parseXML data %) col_names)
       (apply i/conj-cols)
       (i/rename-cols col_replace)))

(def annotation (annotate data col_names replace_col_names))

(defn str>integer [s]
  "returnes the first string number as an integer"
  (-  (Integer. (re-find  #"\d+" s )) 1254))

(defn str>integercol [col_name_in col_name_out file]
  "adds integers columns "
  (->> file
       (i/add-derived-column
        col_name_out
        [col_name_in]
        #(str>integer %))))

(def integrate_select (->> annotation
                    (str>integercol :starts :starts_int)
                    (str>integercol :ends   :ends_int)))

(def integrated (i/$where {:feature {:$eq "mRNA"}} integrate_select))

(def forwerd_seq
  (let [file (i/$where (i/$fn [  starts_int ends_int] 
                              (< starts_int ends_int)) integrated)]
    (->> file
         (i/add-derived-column
          :forwerd_range
          [:starts_int :ends_int]
          #(range %1 %2)))))

(def complement_seq
  (let [file (i/$where (i/$fn [  starts_int ends_int] 
                              (> starts_int ends_int)) integrated)]
    (->> file
         (i/add-derived-column
          :complement_range
          [:ends_int :starts_int]
          #(range %1 %2)))))

(defn bool-vec [vec_size bool]
  "creates a vector of false at vec-size"
  (vec (take vec_size (repeat bool))))

(def hi (vec ( flatten (i/$ :forwerd_range forwerd_seq))))

(defn upd-vec [input-vector ids new-values]
  "interleavs all trues into vector"
  (apply assoc input-vector (interleave ids new-values)))

;(def incar (ii/read-dataset "/home/yosh/datafiles/incanted" :header true))

;(def gi (upd-vec (bool-vec (i/nrow incar) false) (vec hi) (bool-vec (count hi) true)))

;(def added  (i/add-column gi incar))


