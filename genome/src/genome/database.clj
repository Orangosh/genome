(ns genome.database)
(require '[clojure.java.io :as io]
         '[incanter.core :as i]
         '[incanter.datasets :as id]
         '[incanter.io :as ii ]
         '[incanter.charts :as c]
         '[incanter.stats :as st]
         '[clojure.string :as s]
         '[clojure.data.csv :as csv])

;;;;;;;;;;;;;;looks if there in a nil to avoid nullpointer in the d&c_count
;;;;;;;;;;;;;;for a seq


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

(defn unite [col_var col_ref col_name pile_set]
                (->> pile_set
                     (i/add-derived-column
                      col_name
                      [col_var col_ref] +)))

(defn pi [T A G C] 
  (let [cov (+ T A G C)]
    (if (>=  cov 2) 
      (double (/(+ (* T A) (* T G)
                   (* T C) (* A G) 
                   (* A C) (* G C))
                (/ (* cov (- cov 1))
                   2)))
      0)))


(defn SFS [ref T A C G]
  (let [f { "T" T "A" A "C" C "G" G}] 
    (i/sum (filter #(not= (f ref) %) [ T A C G]))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn create-db [file_in file_out]
  (println "Opening a TSV file") 
  
  (def pileup (ii/read-dataset file_in :header false :delim \tab))
  ; "/home/yosh/datafiles/02-519-Pb_CMV_S18whole/output/mpuf"
  (println "Renameing colums")
  
  (def renamed (i/rename-cols
                {:col0     :r_seq
                 :col1     :loc
                 :col2     :ref
                 :col3     :cov
                 :col4     :reads
                 :col5     :qual}
                pileup))
  
  (println "Separating-snp")
  
  (def  seperated-snp (->> renamed
                           (i/add-derived-column
                            :SNPs
                            [:reads]
                            #(sep-snp %))))
  
  (def scrubed
    (->> seperated-snp
         (i/$ [:r_seq :loc :ref :cov :SNPs])))
  
  (println "Adds column of nuc")
  
  (def mapped 
    (->> scrubed
         (i/add-derived-column
          :snap
          [:SNPs]
          #(create-map %1))))

  (def mapped2
    (->> mapped
         (i/$ [:r_seq :loc :ref :cov :snap])))

  (def collumned (->> mapped2
                   (add-col \A) (add-col \a)
                   (add-col \T) (add-col \t)
                   (add-col \C) (add-col \c)
                   (add-col \G) (add-col \g)
                   (add-col \*)))
  
  (def reunited
    (->> (unite \A \a :Aun collumned)
         (unite \T \t :Tun)
         (unite \C \c :Cun)
         (unite \G \g :Gun))) 

  (def finalized
    (->> reunited
         (i/$ [:r_seq :loc :ref :cov :Aun :Tun :Cun :Gun ])))

  (def pied 
    (->> finalized
         (i/add-derived-column
          :pi
          [:Tun :Aun :Gun :Cun]
          #(pi %1 %2 %3 %4))))

  (def SFSd
    (->> pied
         (i/add-derived-column
          :sfs
          [:ref :Tun :Aun :Cun :Gun]
          #(SFS  %1 %2 %3 %4 %5))))
  
  (def nucleotide_diversity (/  (i/sum (i/$ :pi pied)) (i/nrow pied)))
  (def segregation_sites (count (filter #(< 0 %) (i/$ :pi pied))))
  (def coverage (/ (i/sum (i/$ :cov pied)) (i/nrow pied)))
  (println "Nucleotide diversity: " nucleotide_diversity)
  (println "Segregating Sites: " segregation_sites)
  (println "Average coverage: " coverage)

  (with-open [f-out (io/writer file_out)]
    (csv/write-csv f-out [(map name (i/col-names SFSd))])
    (csv/write-csv f-out (i/to-list SFSd))))
