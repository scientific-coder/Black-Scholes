;;; (C) scientific-coder 2011 GPL v3 or later
;;; for now used in incanter console only 6x slower than C++ single core and 10x than multicore !
;; requires incanter/parallelcolt cf project.clj for lein
;;; http://incanter.org/docs/parallelcolt/api/cern/jet/random/tdouble/Normal.html

;;; :DONE: implement formula with cdf-normal
;;; :DONE: multicore
;;; :DONE: move away from incanter :( for 1.3 type hinting :)
(ns scientific-coder.black-scholes
  (:import  java.lang.Math)
  (:import  cern.jet.random.tdouble.Normal)
  (:import  cern.jet.random.tdouble.engine.DoubleMersenneTwister)
  )
(set! *warn-on-reflection* true)
(set! *unchecked-math* true) 
(defn montecarlo-put-call [^double strike-price ^double stock-price]
  (let [n-steps 365
        n-trials 200000
        years-to-maturity 1.
        risk-free-rate 0.03
        volatility 0.2
        coeff (Math/sqrt (/ years-to-maturity n-steps))
        black-scholes (letfn [(cdf-normal
                               [x]
                               (.cdf (Normal. 0. 1. (DoubleMersenneTwister. ) ) x))]
                        (let [d1 (/ (+ (Math/log (/ stock-price strike-price))
                                       (* (+ risk-free-rate (/ (* volatility volatility) 2.))
                                          years-to-maturity))
                                    (* volatility (Math/sqrt years-to-maturity)))
                              d2 (- d1 (* volatility (Math/sqrt years-to-maturity)))
                              call-price (- (* stock-price (cdf-normal d1))
                                            (* strike-price
                                               (Math/exp (* (- risk-free-rate) years-to-maturity))
                                               (cdf-normal d2)))
                              put-price (- (* strike-price
                                              (Math/exp (* (- risk-free-rate) years-to-maturity))
                                              (cdf-normal (- d2)))
                                           (* stock-price (cdf-normal (- d1))))]
                          [put-price call-price]))
        _ (print "formula:" black-scholes)]
    (letfn [(f
             [value]
             (/ (* value (Math/exp (* (- risk-free-rate)
                                      years-to-maturity)))
                n-trials))
            (compute-step
             ^:static ^double [^double price ^double rnd-number]
             (* price (+ 1.
                         (* rnd-number coeff)
                         (* risk-free-rate (/ years-to-maturity n-steps)))))
            
            (make-future-workers
             []
             (let [seed (atom 0)
                   n-iters-todo (atom n-trials)]
               (fn []
                 (future
                   (let [double-normal-rng (Normal. 0. volatility
                                                    (DoubleMersenneTwister. (swap! seed inc)))]
                     (loop [put-value  0.
                            call-value  0.]
                       (if (neg? (swap! n-iters-todo dec))
                         [put-value call-value]
                         (let [delta-price (- strike-price
                                              (reduce compute-step stock-price
                                                      (repeatedly n-steps
                                                                  #(.nextDouble double-normal-rng
                                                                                0. volatility))))
                               ]
                           (if (pos? delta-price)
                             (recur  (+ put-value delta-price) call-value)
                             (recur  put-value (- call-value delta-price)))))))))))]
      (map f (reduce #(map + % (deref %2))
                     [0. 0.]
                     (doall (repeatedly (.. Runtime getRuntime availableProcessors)
                                        (make-future-workers))))))))
(time (montecarlo-put-call 110 100))
