(ns genome.coredraft
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
           [genome.dna2aa     :as da ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MODULES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def none
  [["-a" "--annotate"  "adds annotations"
    :default nil]
   ["-d" "--dna2aa"   ""
    :default false]
   ["-h" "--help"     ""]
   ["-c" "--consvar"  "Deals with variants and concensus calling"
    :default 0]])

(def analyze
  [["-a" "--annotate"  "adds annotations"
    :default nil]
   ["-d" "--dna2aa"   ""
    :default false]
   ["-h" "--help"     ""]
   ["-c" "--consvar"  "Deals with variants and concensus calling"
    :default 0]])

(def database
  [["-a" "--annotate"  "adds annotations"
    :default nil]
   ["-d" "--dna2aa"   ""
    :default false]
   ["-h" "--help"     ""]
   ["-c" "--consvar"  "Deals with variants and concensus calling"
    :default 0]])

(def pops
  [["-a" "--annotate"  "adds annotations"
    :default nil]
   ["-d" "--dna2aa"   ""
    :default false]
   ["-h" "--help"     ""]
   ["-c" "--consvar"  "Deals with "
    :default 0]])

(def stats
  [["-a" "--annotate"  "adds annotations"
    :default nil]
   ["-d" "--dna2aa"   ""
    :default false]
   ["-h" "--help"     ""]
   ["-c" "--consvar"  "Deals with variants and concensus calling"
    :default 0]])

(def view
  [["-a" "--annotate"  "adds annotations"
    :default nil]
   ["-d" "--dna2aa"   ""
    :default false]
   ["-h" "--help"     ""]
   ["-c" "--consvar"  "Deals with variants and concensus calling"
    :default 0]])
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;USAGE MASSAGE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn usage [options-summary first-argument]
  (->> ["Welcome to genome" 
        "This program try to deal with viral sequencing output."
        ""
        "Usage: genome module [options]"
        ""
        "Actions:"
        "  analyse    Analyses different sequences"
        "  database   Print a server's status"
        "  pop        analyses populations"
        "  stats      creates statistics"
        "  view       Visualization instrument"
        ""
        "Options for " first-argument
        options-summary
        ""        
        "Manual page at https://github.com/Orangosh/genome.git"]
       (s/join \newline)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARGUMENTS VALIDATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn error-msg [errors]
  (str "The following errors occurred while parsing your command:\n\n"
       (s/join \newline errors)))

(defn validate-args
  "Validate command line arguments. Either return a map indicating the program
  should exit (with a error message, and optional ok status), or a map
  indicating the action the program should take and the options provided."
  [& args]
  (let [{:keys [options arguments errors summary]}
        (cli/parse-opts args (eval (symbol (first args))))]
    (cond
      (:help options) {:exit-message (usage summary (first arguments)) :ok? true}
      errors          {:exit-message (error-msg errors)}
      (and (= 1 (count arguments))
           (#{"none" "analyze" "database" "pops" "stats" "wiew"} (first arguments)))
                      {:action (first arguments) :options options}
      :else           {:exit-message (usage summary (first arguments))})))

(defn exit [status msg]
  (println msg)
  (System/exit status))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAIN FUNCTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn -main [module file_in file_out & args]
  (let [{:keys [action options exit-message ok?]} (validate-args args)]
    (if exit-message
      (exit (if ok? 0 1) exit-message)
      (case action   ;first argument
        "analyse"                           ;Analyses different sequences
        "database"   (gd/create-db options) ;Creates a database
        "pops"                              ;analyses populations
        "stats"                             ;creates statistics
        "view"                              ;Visualization instrument
        ))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TESTING PIPELINE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn -create_db [file_in file_out]
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

