(ns genome.database
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FUNCTIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn drop-parse [st]
;st start with number- and drops number and strign folowing in number size
  (let [s (Integer. (re-find  #"\d+" st ))] 
    (apply str (drop (+ (count (str s)) s) st))))

(defn sep-snp [sq]
 (if (not= nil sq)
   (let [seq (-> sq
                 (s/replace #"\^." "") ;droping problematic values before parsing
                 (s/replace #"\-" "")
                 (s/replace #"\+" ""))]
     (->> seq 
          (apply str)
          (re-seq  #"\d+[^\d]*") ;seq of of string that start with numbers
          (map drop-parse)
          (apply str (apply str (re-seq #"^[^\d]+" seq))))) ;adds the first string
   "!"))

(defn create-map [seq]
  (let [sym {\A 0 \a 0 \T 0 \t 0
             \C 0 \c 0 \G 0 \g 0
             \. 0 \, 0 \* 0}
        freq (frequencies seq)]
    (merge sym freq)))

(defn merge-all [col_name ref merged]
  (cond
    (= (str col_name) ref)
      (merged \.)
      (= (str col_name) (s/lower-case ref))
      (merged \,)
      :else (merged col_name)))

(defn add-col [col_name dbase]
  (->> dbase
       (i/add-derived-column
        col_name
        [:ref :snap]
         #(merge-all col_name %1 %2))))

(defn unite [col_fwd col_bwd col_name strand_filter pile_set]
  "stand_filter- number 0-1- ignores variants that are (90% default) dominant
0.0 means no filter"
                (->> pile_set
                     (i/add-derived-column
                      col_name
                      [col_fwd col_bwd]
                      #(if (and (< 0.0 %1) (< 0.0 %2))
                         (if (and (> (double (/ %1 %2)) strand_filter)
                                  (> (double (/ %2 %1)) strand_filter))
                           (+ %1 %2)
                           0)
                         0))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;OPPERATION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn create-db [file & {:keys [strand_filter]
                         :or {strand_filter 0.0}}]
  (let [nucleotides [\A \a \T \t \C \c \G \g \*]]
    (->> (ii/read-dataset file :header false :delim \tab)
         (i/rename-cols
          {:col0 :r_seq
           :col1 :loc
           :col2 :ref
           :col3 :cov
           :col4 :reads
           :col5 :qual})
         (i/add-derived-column :SNPs
          [:reads]
          #(sep-snp %))
         (i/add-derived-column :snap
          [:SNPs]
          #(create-map %))
         (add-col \A) (add-col \a)
         (add-col \T) (add-col \t)
         (add-col \C) (add-col \c)
         (add-col \G) (add-col \g)
         (add-col \*)
         (unite \A \a :Aun strand_filter)
         (unite \T \t :Tun strand_filter)
         (unite \C \c :Cun strand_filter)
         (unite \G \g :Gun strand_filter)
         (i/$ (vec (flatten [:r_seq :loc :ref  :cov
                             :Aun  :Tun  :Cun  :Gun
                             nucleotides]))))))

