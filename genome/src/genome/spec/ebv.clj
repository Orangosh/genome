(ns genome.spec.ebv
  (require [clojure.java.io   :as io ]
           [incanter.core     :as i  ]
           [incanter.io       :as ii ]
           [incanter.stats    :as st ]
           [incanter.charts   :as c  ]
           [genome.dna2aa     :as da ]
           [genome.stats      :as gs ]
           [genome.consvar    :as cv ]
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE FILE LOCATIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn ebv-loc []
  (def L44-Pb  (str home "ebv/incanted_files/344-Pb.inc" ))   
  (def L44-S2a (str home "ebv/incanted_files/344-S2a.inc"))  
  (def L44-S2b (str home "ebv/incanted_files/344-S2b.inc"))  
  (def L44-S2c (str home "ebv/incanted_files/344-S2c.inc"))   
  (def L44-S3  (str home "ebv/incanted_files/344-S3.inc" ))
  
  (def L25-Ma  (str home "ebv/incanted_files/525-Ma.inc" ))  
  (def L25-Mc  (str home "ebv/incanted_files/525-Mc.inc" ))  
  (def L25-Pa  (str home "ebv/incanted_files/525-Pa.inc" ))  
  (def L25-S1a (str home "ebv/incanted_files/525-S1a.inc"))  
  (def L25-S1b (str home "ebv/incanted_files/525-S1b.inc"))  
  (def L25-S1c (str home "ebv/incanted_files/525-S1c.inc"))  
  (def L25-S1d (str home "ebv/incanted_files/525-S1d.inc"))  
  (def L25-S2a (str home "ebv/incanted_files/525-S2a.inc"))  
  (def L25-S2b (str home "ebv/incanted_files/525-S2b.inc"))  
  (def L25-S3a (str home "ebv/incanted_files/525-S3a.inc"))  
  (def L25-S3b (str home "ebv/incanted_files/525-S3b.inc"))  
  (def L25-S3c (str home "ebv/incanted_files/525-S3c.inc"))  

  (def L38-Ma  (str home "ebv/incanted_files/538-Ma.inc" ))   
  (def L38-Pb  (str home "ebv/incanted_files/538-Pb.inc" ))   
  (def L38-S1a (str home "ebv/incanted_files/538-S1a.inc"))  

  (def L40-Ma  (str home "ebv/incanted_files/540-Ma.inc" ))  
  (def L40-Mb  (str home "ebv/incanted_files/540-Mb.inc" ))  
  (def L40-Mc  (str home "ebv/incanted_files/540-Mc.inc" ))  
  (def L40-Pa  (str home "ebv/incanted_files/540-Pa.inc" ))
  (def L40-Pc  (str home "ebv/incanted_files/540-Pc.inc" ))
  (def L40-S1a (str home "ebv/incanted_files/540-S1a.inc"))
  (def L40-S1b (str home "ebv/incanted_files/540-S1b.inc")))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET SET
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-set [file cov]
  "open an csv.inc file"
  (->> (ii/read-dataset file :header true)
       (i/$where (i/$fn [depth] (< cov depth)))))
(def m-get-set (memoize get-set))

(defn ebv-sets []
  (ebv-loc)
  (def samples {:E44-Pb  (m-get-set L44-Pb  0)   
                :E44-S2a (m-get-set L44-S2a 0)  
                :E44-S2b (m-get-set L44-S2b 0)  
                :E44-S2c (m-get-set L44-S2c 0)   
                :E44-S3  (m-get-set L44-S3  0)
                
                :E25-Ma  (m-get-set L25-Ma  0)  
                :E25-Mc  (m-get-set L25-Mc  0)  
                :E25-Pa  (m-get-set L25-Pa  0)  
                :E25-S1a (m-get-set L25-S1a 0)  
                :E25-S1b (m-get-set L25-S1b 0)  
                :E25-S1c (m-get-set L25-S1c 0)  
                :E25-S1d (m-get-set L25-S1d 0)  
                :E25-S2a (m-get-set L25-S2a 0)  
                :E25-S2b (m-get-set L25-S2b 0)  
                :E25-S3a (m-get-set L25-S3a 0)  
                :E25-S3b (m-get-set L25-S3b 0)  
                :E25-S3c (m-get-set L25-S3c 0)  

                :E38-Ma  (m-get-set L38-Ma  0)   
                :E38-Pb  (m-get-set L38-Pb  0)   
                :E38-S1a (m-get-set L38-S1a 0)  

                :E40-Ma  (m-get-set L40-Ma  0)  
                :E40-Mb  (m-get-set L40-Mb  0)  
                :E40-Mc  (m-get-set L40-Mc  0)  
                :E40-Pa  (m-get-set L40-Pa  0)
                :E40-Pc  (m-get-set L40-Pc  0)
                :E40-S1a (m-get-set L40-S1a 0)
                :E40-S1b (m-get-set L40-S1b 0)}))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;FOR PCA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn le-filter [file & {:keys [m d]
                          :or   {m 0.05
                                 d 35.0}}]
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
                   depth17 depth18 depth19
                   depth20 depth21
                   depth22 depth23 depth24
                   depth25 depth26 depth27
                   ] 
                  (and (not= nil depth1)
                       (not= nil depth2)  (not= nil depth3)  (not= nil depth4)
                       (not= nil depth5)  (not= nil depth6)  (not= nil depth7)
                       (not= nil depth8)  (not= nil depth9)
                       (not= nil depth10) (not= nil depth11)
                       (not= nil depth12)
                       (not= nil depth13) (not= nil depth14)
                       (not= nil depth15) (not= nil depth16)
                       (not= nil depth17) (not= nil depth18) (not= nil depth19)
                       (not= nil depth20) (not= nil depth21)
                       (not= nil depth22) (not= nil depth23) (not= nil depth24)
                       (not= nil depth25) (not= nil depth26) (not= nil depth27)
                       ))) 
       (i/$where (i/$fn
                  [minfr1
                   minfr2  minfr3  minfr4
                   minfr5  minfr6  minfr7
                   minfr8  minfr9
                   minfr10 minfr11
                   minfr12
                   minfr13 minfr14
                   minfr15 minfr16
                   minfr17 minfr18 minfr19
                   minfr20 minfr21
                   minfr22 minfr23 minfr24
                   minfr25 minfr26 minfr27
                   ]
                  (or  (> minfr1 m)
                       (> minfr2 m)  (> minfr3 m)  (> minfr4 m)
                       (> minfr5 m)  (> minfr6 m)  (> minfr7 m)
                       (> minfr8 m)  (> minfr9 m)
                       (> minfr10 m) (> minfr11 m)
                       (> minfr12 m)
                       (> minfr13 m) (> minfr14 m)
                       (> minfr15 m) (> minfr16 m)
                       (> minfr17 m) (> minfr18 m) (> minfr19 m)
                       (> minfr20 m) (> minfr21 m)
                       (> minfr22 m) (> minfr23 m) (> minfr24 m)
                       (> minfr25 m) (> minfr26 m) (> minfr27 m)
                       ))) 
       (i/$where (i/$fn
                  [depth1
                   depth2  depth3  depth4
                   depth5  depth6  depth7
                   depth8  depth9
                   depth10 depth11
                   depth12
                   depth13 depth14
                   depth15 depth16
                   depth17 depth18 depth19
                   depth20 depth21
                   depth22 depth23 depth24
                   depth25 depth26  depth27
                   ]
                  (and (> depth1  d)
                       (> depth2  d) (> depth3  d) (> depth4  d)
                       (> depth5  d) (> depth6  d) (> depth7  d)
                       (> depth8  d) (> depth9  d)
                       (> depth10 d) (> depth11 d)
                       (> depth12 d)
                       (> depth13 d) (> depth14 d)
                       (> depth15 d) (> depth16 d)
                       (> depth17 d) (> depth18 d) (> depth19 d)
                       (> depth20 d) (> depth21 d)
                       (> depth22 d) (> depth23 d) (> depth24 d)
                       (> depth25 d) (> depth26 d) (> depth27 d)
                       )))
       (i/$ [:minfr1
             :minfr2  :minfr3  :minfr4
             :minfr5  :minfr6  :minfr7
             :minfr8  :minfr9
             :minfr10 :minfr11
             :minfr12
             :minfr13 :minfr14
             :minfr15 :minfr16
             :minfr17 :minfr18 :minfr19
             :minfr20 :minfr21
             :minfr22 :minfr23 :minfr24
             :minfr25 :minfr26 :minfr27
             ])))

#_(genome.dreduction/save-mat (le-filter mat) "/mnt/data/ebv/ebv-matrix")

(defn EBV-SVD-depth [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 7 8 12 18 20 23] 0 projection)
                        (i/$ [0 7 8 12 18 20 23] 1 projection)
                        :title "EBV by depth"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [1 5 6 10 15 19 24 25] 0 projection)
                      (i/$ [1 5 6 10 15 19 24 25] 1 projection))
        (c/add-points (i/$ [2 3 4 9 11 13, 14 16 17 21 22 26] 0 projection)
                      (i/$ [2 3 4 9 11 13, 14 16 17 21 22 26] 1 projection))
        (i/view))))


(defn EBV-SVD-family [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ (range 0 5) 0 projection)
                        (i/$ (range 0 5) 1 projection)
                        :title "EBV by family"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ (range 5 17)  0 projection)
                      (i/$ (range 5 17)  1 projection))
        (c/add-points (i/$ (range 17 20) 0 projection)
                      (i/$ (range 17 20) 1 projection))
        (c/add-points (i/$ (range 20 26) 0 projection)
                      (i/$ (range 20 26) 1 projection))
        (i/view))))

(defn EBV-SVD-primary [file]
  (let [projection (->> (ii/read-dataset file :header false) i/to-matrix)]
    (-> (c/scatter-plot (i/$ [0 7 8 12 17 19 22] 0 projection)
                        (i/$ [0 7 8 12 17 19 23] 1 projection)
                        :title "EBV primary vs mothers and siblings"
                        :x-label "Dimension 1"
                        :y-label "Dimension 2")
        (c/add-points (i/$ [2 3 4 9 11 13 15 16 20 21 24 25] 0 projection)
                      (i/$ [2 3 4 9 11 13 15 16 20 21 24 25] 1 projection))
        (c/add-points (i/$ [1 5 6 10 14 18 23 24] 0 projection)
                      (i/$ [1 5 6 10 14 18 23] 1 projection))
        (i/view))))




