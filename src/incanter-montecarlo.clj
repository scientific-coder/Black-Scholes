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
                        [put-price call-price])]
    (letfn [(f
             [value]
             (/ (* value (java.lang.Math/exp (* (- risk-free-rate)
                                                years-to-maturity)))
                num-trials))
            (compute-step
             [rnd-nums-init]
             (loop [rnd-nums rnd-nums-init
                    i n-steps
                    result stock-price]
               (if (== i 0)
                 [result rnd-nums]
                 (recur (next rnd-nums)
                        (dec i)
                        (* result (+ 1. (* (first rnd-nums)
                                           (java.lang.Math/sqrt (/ years-to-maturity
                                                                   n-steps)))
                                     (* risk-free-rate
                                        (/ years-to-maturity n-steps))))))))]
      (loop [trial num-trials
             put-value 0.
             call-value 0.
             random-nbs (sample-normal (* n-steps num-trials) :mean 0. :sd volatility)]
        (if (== trial 0)
          (do (print "formula gave:" black-scholes) (map f [put-value call-value]))
          (let [[computed-price next-random-nbs] (compute-step random-nbs)
                delta-price (- strike-price computed-price)
                [next-put next-call] (if (pos? delta-price)
                                       [(+ put-value delta-price) call-value]
                                       [put-value (- call-value delta-price)])]
            ( recur (dec trial) next-put next-call next-random-nbs)))))))
(time (montecarlo-put-call 110 100))
