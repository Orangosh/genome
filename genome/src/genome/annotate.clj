(ns genome.annotate
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.stats :as st]
           [clojure.xml :as xml]
           [clojure.zip :as zip]))

;; wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_006273.2&retmode=xml"
;; mv efetch.fcgi?db=nuccore\&id=NC_006273.2\&retmode=xml merlin.xm

(def data "/home/yosh/datafiles/merlin.xml")

(def col_names {:GBFeature_key :feature
                :GBInterval_from :starts
                :GBInterval_to :ends})

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
  (->> (map #(parseXML data %) (vec (keys col_names)))
       (apply i/conj-cols)
       (i/rename-cols col_replace)))

(def annotation (annotate data col_names replace_col_names))

(defn str>integer [s]
  "returnes the first string number as an integer"
  (Integer. (re-find  #"\d+" s )))

(defn str>integercol [col_name_in col_name_out file]
  ""
  (->> file
       (i/add-derived-column
        col_name_out
        [col_name_in]
        #(str>integer %))))

(def integered (->> annotation
                    (str>integercol :starts :starts_int)
                    (str>integercol :ends   :ends_int)))

(def forwerd_seq
  (let [file (i/$where (i/$fn [  starts_double ends_double] 
                              (< starts_double ends_double)) i)]
    (->> file
         (i/add-derived-column
          :forwerd_range
          [:starts_int :ends_int]
          #(range %1 %2)))))

(def complement_seq
  (let [file (i/$where (i/$fn [  starts_double ends_double] 
                              (> starts_double ends_double)) i)]
    (->> file
         (i/add-derived-column
          :complement_range
          [:ends_int :starts_int]
          #(range %1 %2)))))

(defn bool-vec [file bool]
  "creates a vector of false at vec-size"
  (let [vec-size (i/nrow file)]
       (vec (take vec-size (repeat bool)))))

(def hi (vec ( flatten (i/$ :complement_range complement_seq))))

(defn upd-vec [input-vector ids new-values]
  ""
  (apply assoc input-vector (interleave ids new-values)))

(def gi (upd-vec (bool-vec 235000 false) hi (bool-vec (count hi) true)))


