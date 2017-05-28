(ns genome.core
  (:gen-class)
  (require [clojure.java.io   :as io ]
           [clojure.string    :as s  ]
           [clojure.data.csv  :as csv]
           [clojure.tools.cli :as cli]
           [incanter.core     :as i  ]
           [incanter.datasets :as id ]
           [incanter.io       :as ii ]
           [incanter.charts   :as c  ]
           [incanter.stats    :as st ]
           [genome.database   :as gd ]
           [genome.stats      :as gs ]
           [genome.pop        :as p  ]
           [genome.consvar    :as cv ]
           [genome.refloc     :as rl ]
           [genome.dna2aa     :as da ]
           [genome.annotate   :as ga ]
           [genome.view       :as v  ]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TESTING PIPELINE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn -main [file_in file_out refset gff3]
  (println "Welcome to clojure- starting incanter")
  (gd/create-db file_in)

  (println "Adding ref location")
  (def refered     (rl/refer-ann refset gd/finalized))

  (println "Adding annotations")
  (def annotated   (ga/annotate gff3 refered))

  (println "First dept adjustment")
  (def un_coved    (cv/calc-coved cv/cov_un annotated))

  (println "Creating first consensus sequence")
  (def fstsensus    (cv/major-allele cv/maj_un un_coved))
  
  (println "Correcting read errors")
  (def minored
  "minor allele-> cv/minorallize/poisson dist-> cv/poissonize"
    (cv/get-minor-allele "minor allele" 0.05 fstsensus))
  
  (println "Corrected dept adjustment")
  (def coved  (cv/calc-coved cv/depth minored))

  (println "Creating second consensus sequence")
  (def majored     (cv/major-allele cv/maj coved))

  (println "Adding consensus negative strand")
  (def neg_majored (da/pos>neg :maj+ :maj- majored))
  
  (println "Creating minor allele sequence")
  (def minored     (cv/minor-allele cv/min_p neg_majored))  

  (println "Adding negative strand")
  (def neg_minored (da/pos>neg :min+ :min- minored))
  
  (println "adding consensus amino acids")
  (def maj_aa      (da/nuc>aa  :maj+ :maj- neg_minored))

  (println "adding consensus amino acids")
  (def min_aa      (da/nuc>aa :min+ :min- maj_aa))  

  (def scrubed2 (i/$ [:r_seq
                      :merlin :ref-loc
                      :gene+ :gene- :CDS+ :CDS- :exon- :exon+
                      :ref :loc :maj_un+
                      :cov :cov_un :depth
                      :maj+ :min+ :maj_aa+ :min_aa+
                      :maj- :min- :maj_aa- :min_aa-
                      :T :A :C :G] min_aa))
  
  (println "Calculating nucleotide diversity")
  (def pied (p/pise p/pi scrubed2))

  (println "Calculating folded allele frequency spectra")
  (def sfsd (p/SFS p/folded-SFS pied))

  (println "Adding orfs")
  (def orfd (da/add-orf-aa sfsd))

  (with-open [f-out (io/writer file_out)]
    (csv/write-csv f-out [(map name (i/col-names orfd))])
    (csv/write-csv f-out (i/to-list orfd))))

(defn ready []
  (-main "/home/yosh/datafiles/experiment/mpileup" 
         "/home/yosh/datafiles/experiment/incanter"
         "/home/yosh/datafiles/Consensus/CMVconsensus/refset.inc"
         "/home/yosh/datafiles/genes/merlin.gff3"))


(defn show[from length]
  (def incanted (ii/read-dataset "/home/yosh/datafiles/incanter" :header true))
  (i/$ (range from (+ from length)) :all incanted))

