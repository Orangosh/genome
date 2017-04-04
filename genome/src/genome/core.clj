(ns genome.core
  (:gen-class)
  (require [clojure.java.io :as io]
           [incanter.core :as i]
           [incanter.datasets :as id]
           [incanter.io :as ii ]
           [incanter.charts :as c]
           [incanter.stats :as st]
           [clojure.string :as s]
           [clojure.data.csv :as csv]
           [genome.database :as gd]
           [genome.stats :as gs]
           [genome.pop :as p]
           [genome.consvar :as cv]
           [genome.dna2aa :as da]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TESTING PIPELINE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn -main [file_in file_out]
  (println "Welcome to clojure- starting incanter")
  (gd/create-db file_in)

  (println "Creating first consensus")
  (def conded (cv/consensus gd/finalized cv/consus_un))

  (println "Correcting read errors")
  (def pois (cv/poissonize 0.05 conded))

  (def scrubed2 (i/$ [:r_seq :loc :ref :consus_un :cov :c_cov
                      :Tpois :Apois :Cpois :Gpois] pois))
  (println "Calculating nucleotide diversity")
  (def pied (p/pise scrubed2 p/pi_pois))

  (println "Calculating folded allele frequency spectra")
  (def sfsd (p/SFS pied p/folded-SFS))

  (println "adding negative strand")
  (def neg_stranded (da/pos>neg sfsd :consus_un :negsus_un))

  (println "adding consensus amino acids")
  (def aaadded (da/nuc>aa neg_stranded :consus_un :negsus_un))

  (println "removing INDELs")
  (def row_cleaned (gs/row-clean aaadded :ref "-"))
  
  (println "SUMMARY STATISTICS:")
  (gs/stat-report sfsd)

  (with-open [f-out (io/writer file_out)]
    (csv/write-csv f-out [(map name (i/col-names aaadded))])
    (csv/write-csv f-out (i/to-list aaadded))))


(defn ready []
  (-main "/home/yosh/datafiles/mpileup" "/home/yosh/datafiles/incanter"))


(defn show[from length]
  (def incanted (ii/read-dataset "/home/yosh/datafiles/incanter" :header true))
  (i/$ (range from (+ from length)) :all incanted))
