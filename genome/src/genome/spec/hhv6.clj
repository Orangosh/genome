(ns genome.spec.hhv6
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [incanter.charts   :as c  ]
           [genome.dna2aa     :as da ]
           [genome.stats      :as gs ]

           [genome.pop        :as p  ]
           [clojure.data.csv  :as csv]
           [genome.view       :as v  ]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;
;; For PC
;;(def home "/home/yosh/data")

;;;;;;;;;;;;;;;;;;;;;;;
;; For Server
(def home "/mnt/data/")


(defn hhv6-loc []
  ;;hhv6-sets 
  (def L37-Pa  (str home "hhv6/incanted_files/537-Pa.inc" ))
  (def L37-S2b (str home "hhv6/incanted_files/537-S2b.inc"))        
  (def L37-Pb  (str home "hhv6/incanted_files/537-Pb.inc" ))
  (def L37-S3  (str home "hhv6/incanted_files/537-S3.inc" ))       
  (def L37-Pc  (str home "hhv6/incanted_files/537-Pc.inc" ))            
  (def L37-S2a (str home "hhv6/incanted_files/537-S2a.inc"))

  (def L42-Mb  (str home "hhv6/incanted_files/542-Mb.inc" ))
  (def L42-Ma  (str home "hhv6/incanted_files/542-Ma.inc" ))
  (def L42-Pa  (str home "hhv6/incanted_files/542-Pa.inc" ))
  (def L42-Pb  (str home "hhv6/incanted_files/542-Pb.inc" ))
  (def L42-Pc  (str home "hhv6/incanted_files/542-Pc.inc" ))
  (def L42-S1a (str home "hhv6/incanted_files/542-S1a.inc"))
  (def L42-S1b (str home "hhv6/incanted_files/542-S1b.inc"))
  (def L42-S1c (str home "hhv6/incanted_files/542-S1c.inc"))

  (def L43-Pb  (str home "hhv6/incanted_files/543-Pb.inc" ))
  (def L43-S1a (str home "hhv6/incanted_files/543-S1a.inc"))  
  (def L43-S1b (str home "hhv6/incanted_files/543-S1b.inc"))

  (def L72-M   (str home "hhv6/incanted_files/572-M.inc"  ))
  (def L72-Pa  (str home "hhv6/incanted_files/572-Pa.inc" )))


(defn get-set [file cov]

  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))


(defn hhv6-sets []
  (hhv6-loc)
  (def H37-Pa  (m-get-set L37-Pa  0))
  (def H37-S2b (m-get-set L37-S2b 0))        
  (def H37-Pb  (m-get-set L37-Pb  0))
  (def H37-S3  (m-get-set L37-S3  0))       
  (def H37-Pc  (m-get-set L37-Pc  0))            
  (def H37-S2a (m-get-set L37-S2a 0))

  (def H42-Mb  (m-get-set L42-Mb  0))
  (def H42-Ma  (m-get-set L42-Ma  0))
  (def H42-Pa  (m-get-set L42-Pa  0))
  (def H42-Pb  (m-get-set L42-Pb  0))
  (def H42-Pc  (m-get-set L42-Pc  0))
  (def H42-S1a (m-get-set L42-S1a 0))
  (def H42-S1b (m-get-set L42-S1b 0))
  (def H42-S1c (m-get-set L42-S1c 0))

  (def H43-Pb  (m-get-set L43-Pb  0))
  (def H43-S1a (m-get-set L43-S1a 0))  
  (def H43-S1b (m-get-set L43-S1b 0))

  (def H72-M   (m-get-set L72-M   0))
  (def H72-Pa  (m-get-set L72-Pa  0)))



(def samples [H37-Pa  H37-Pb  H37-Pc
                            H37-S2a H37-S2b H37-S3
                            H42-Ma  H42-Mb
                            H42-Pa  H42-Pb  H42-Pc
                            H42-S1a H42-S1b H42-S1c
                            H43-Pb  H43-S1a H43-S1b
              H72-M   H72-Pa])



(defn le-filter [file & {:keys [m d]
                          :or   {m 0.01
                                 d 13.0}}]
  (->> file
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  depth4
                   depth5  depth6  depth7
                   depth8  depth9
                   depth10 depth11
                   depth12 depth13
                   depth14  
                   depth15 depth16
                   depth17 depth18 depth19] 
                  (and (not= nil depth1)
                       (not= nil depth2)  (not= nil depth3)  (not= nil depth4)
                       (not= nil depth5)  (not= nil depth6)  (not= nil depth7)
                       (not= nil depth8)  (not= nil depth9)
                       (not= nil depth10) (not= nil depth11)
                       (not= nil depth12)
                       (not= nil depth13) (not= nil depth14)
                       (not= nil depth15) (not= nil depth16)
                       (not= nil depth17) (not= nil depth18) (not= nil depth19)))) 
       (i/$where (i/$fn
                  [minfr1
                   minfr2  minfr3  minfr4
                   minfr5  minfr6  minfr7
                   minfr8  minfr9
                   minfr10 minfr11
                   minfr12
                   minfr13 minfr14
                   minfr15 minfr16
                   minfr17 minfr18 minfr19]
                  (or  (> minfr1 m)
                       (> minfr2 m)  (> minfr3 m)  (> minfr4 m)
                       (> minfr5 m)  (> minfr6 m)  (> minfr7 m)
                       (> minfr8 m)  (> minfr9 m)
                       (> minfr10 m) (> minfr11 m)
                       (> minfr12 m)
                       (> minfr13 m) (> minfr14 m)
                       (> minfr15 m) (> minfr16 m)
                       (> minfr17 m) (> minfr18 m) (> minfr19 m)))) 
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  depth4
                   depth5  depth6  depth7
                   depth8  depth9
                   depth10 depth11
                   depth12
                   depth13 depth14
                   depth15 depth16
                   depth17 depth18 depth19]
                  (and (> depth1  d)
                       (> depth2  d) (> depth3  d) (> depth4  d)
                       (> depth5  d) (> depth6  d) (> depth7  d)
                       (> depth8  d) (> depth9  d)
                       (> depth10 d) (> depth11 d)
                       (> depth12 d)
                       (> depth13 d) (> depth14 d)
                       (> depth15 d) (> depth16 d)
                       (> depth17 d) (> depth18 d) (> depth19 d))))
       (i/$ [:minfr1
             :minfr2  :minfr3  :minfr4
             :minfr5  :minfr6  :minfr7
             :minfr8  :minfr9
             :minfr10 :minfr11
             :minfr12
             :minfr13 :minfr14
             :minfr15 :minfr16
             :minfr17 :minfr18 :minfr19])))

(defn HHV6-SVD-primary [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 1 2 8 9 10 14 18] 0 projection)
                        (i/$ [0 1 2 8 9 10 14 18] 1 projection)
                        :title "HHV6 primary vs mother vs sibling"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [3 4 5 11 12 13 15 16]  0 projection)
                      (i/$ [3 4 5 11 12 13 15 16]  1 projection))
        (c/add-points (i/$ [6 7 17] 0 projection)
                      (i/$ [6 7 17] 1 projection))
        (i/view))))


(defn HHV6-SVD-family [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 1 2 3 4 5] 0 projection)
                        (i/$ [0 1 2 3 4 5] 1 projection)
                        :title "hhv6 family"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [6 7 8 9 10 11 12 13]  0 projection)
                      (i/$ [6 7 8 9 10 11 12 13]  1 projection))
        (c/add-points (i/$ [14 15 16] 0 projection)
                      (i/$ [14 15 16] 1 projection))
        (c/add-points (i/$ [17 18] 0 projection)
                      (i/$ [17 18] 1 projection))
        (i/view))))

(defn HHV6-SVD-depth [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 1 2 9 10 14 15 16] 0 projection)
                        (i/$ [0 1 2 9 10 14 15 16] 1 projection)
                        :title "hhv6 depth"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [3 4 5 6 7 8 11 12 17 18]  0 projection)
                      (i/$ [3 4 5 6 7 8 11 12 17 18]  1 projection))
        (i/view))))
