# for now used in incanter console only 6x slower than C++ single core and 18x than multicore !
(use '(incanter core stats))

(comment 
  (cern.jet.random.tdouble.Normal/staticNextDouble mean sd)
  (new cern.jet.random.tdouble.Normal 10 5 (new cern.jet.random.tdouble.engine.DoubleMersenneTwister )))
#http://incanter.org/docs/parallelcolt/api/cern/jet/random/tdouble/Normal.html

# :DONE: implement formula with cdf-normal

(defn montecarlo-put-call [strike-price stock-price]
  (let [n-steps 365
        num-trials 200000
        years-to-maturity 1.
        risk-free-rate 0.03
        volatility 0.2
        black-scholes (let [d1 (/ (+ (java.lang.Math/log (/ stock-price strike-price))
                                     (* (+ risk-free-rate (/ (* volatility volatility) 2.))
                                        years-to-maturity))
                                  (* volatility (java.lang.Math/sqrt years-to-maturity)))
                            d2 (- d1 (* volatility (java.lang.Math/sqrt years-to-maturity)))
                            call-price (- (* stock-price (cdf-normal d1))
                                          (* strike-price (java.lang.Math/exp (* (- risk-free-rate) years-to-maturity))
                                             (cdf-normal d2)))
                            put-price (- (* strike-price
                                            (java.lang.Math/exp (* (- risk-free-rate) years-to-maturity))
                                            (cdf-normal (- d2)))
                                         (* stock-price (cdf-normal (- d1))))]
                        [put-price call-price])
        _ (print "formula:" black-scholes)]
    (letfn [(f
             [value]
             (/ (* value (java.lang.Math/exp (* (- risk-free-rate)
                                                years-to-maturity)))
                num-trials))
            (compute-step
             [price rnd-number]
             (* price (+ 1.
                         (* rnd-number (java.lang.Math/sqrt (/ years-to-maturity n-steps)))
                         (* risk-free-rate (/ years-to-maturity n-steps)))))
            (combine
             [[put-value call-value] rnd-nums]
             (let [delta-price (- strike-price (reduce compute-step stock-price rnd-nums))]
               (if (pos? delta-price)
                 [(+ put-value delta-price) call-value]
                 [put-value (- call-value delta-price)])))]
      (map f (reduce combine [0. 0.]
                     (partition n-steps (sample-normal (* n-steps num-trials)
                                                       :mean 0. :sd volatility)))))))
(time (montecarlo-put-call 110 100))
