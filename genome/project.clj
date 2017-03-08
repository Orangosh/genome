(defproject genome "0.1.1-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [clojure-csv "2.0.2"]
                 [incanter "1.5.7"]     
                 [org.clojure/data.csv "0.1.3"]
                 [semantic-csv "0.1.0"]
                 [yieldbot/vizard "0.1.0"]
                 [cljam "0.1.5"]]
  :jvm-opts ["-Xmx15g"]
  :main genome.core
  :aot [genome.core])
