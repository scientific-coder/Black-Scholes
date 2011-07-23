;;; (C) scientific-coder 2011 GPL v3 or later
;;; for now used in incanter console only 6x slower than C++ single core and 18x than multicore !
;;
;;; http://incanter.org/docs/parallelcolt/api/cern/jet/random/tdouble/Normal.html

;;; :DONE: implement formula with cdf-normal
;;; :DONE: multicore
;;; :DONE: move away from incanter :( for 1.3 type hinting :)
(ns scientific-coder.black-scholes
  (:import  java.lang.Math)
  (:import  cern.jet.random.tdouble.Normal)
  (:import  cern.jet.random.tdouble.engine.DoubleMersenneTwister)
  )
(defn montecarlo-put-call [strike-price stock-price]
  (let [n-steps 365
        n-trials 200
        years-to-maturity 1.
        risk-free-rate 0.03
        volatility 0.2
        black-scholes (letfn [(cdf-normal
                               [x]
                               (.cdf (Normal. 0. 1. (DoubleMersenneTwister. ) ) x))]
                        (let [d1 (/ (+ (Math/log (/ stock-price strike-price))
                                       (* (+ risk-free-rate (/ (* volatility volatility) 2.))
                                          years-to-maturity))
                                    (* volatility (Math/sqrt years-to-maturity)))
                              d2 (- d1 (* volatility (Math/sqrt years-to-maturity)))
                              call-price (- (* stock-price (cdf-normal d1))
                                            (* strike-price (Math/exp (* (- risk-free-rate) years-to-maturity))
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
             [price rnd-number]
             (* price (+ 1.
                         (* rnd-number (Math/sqrt (/ years-to-maturity n-steps)))
                         (* risk-free-rate (/ years-to-maturity n-steps)))))
            
            (make-future-workers
             []
             (let [seed (atom 0)
                   n-iters-todo (atom n-trials)]
               (fn []
                 (future
                   (let [double-normal-rng (Normal. 0. volatility
                                                (DoubleMersenneTwister. (swap! seed inc)))]
                     (loop [put-value 0.
                            call-value 0.]
                       (if (neg? (swap! n-iters-todo dec))
                         [put-value call-value]
                         (let [delta-price (- strike-price (reduce compute-step stock-price
                                                                   (repeatedly n-steps #(.nextDouble double-normal-rng 0. volatility))))
                               [next-put-value next-call-value] (if (pos? delta-price)
                                                                  [(+ put-value delta-price) call-value]
                                                                  [put-value (- call-value delta-price)])]
                           (recur next-put-value next-call-value)))))))))]
      (map f (reduce #(map + % (deref %2))
                     [0. 0.]
                     (doall (repeatedly (.. Runtime getRuntime availableProcessors)
                                        (make-future-workers))))))))
(time (montecarlo-put-call 110 100))
