// Copyright Bernard Hugueney 2011, Licence GPL v3
#include <atomic>
#include <future>
#include <random>
#include <tr1/random> // needed for std::tr1::variate_generator
#include <numeric>
#include <boost/progress.hpp>
#include <boost/math/distributions/normal.hpp>
#include <iostream>
#include <time.h>
#include <chrono>

/* Work splitting between thread easier because we work from rng output instead of a collection of data.
   So we do not have to worry about memory locality. But we do have to worry about the rng state : we do not want to share
   it between threads (to avoid concurrency issues), but we need to ensure that each thread has as unique rng state 
   : we cannot copy it from one thread to another because each would then work on the same pseudo random output which would defeat
   the whole purpose of having multiple results !
   We resolve this by sharing a seed that each thread will use upon creation to get a unique state.
 */
//g++ cxx0x-montecarlo.cxx -o cxx0x-montecarlo -std=c++0x -Wall -O4 -march=native -lpthread
struct seeder_t{
  unsigned int operator()() const
  { return static_cast<unsigned int>(time(0))+ (++c); }
  static std::atomic<unsigned int> c; // current code does not need atomic<>
};
std::atomic<unsigned int> seeder_t::c(0);

// each thread will get it's own private engine, and each each will have to be re seeded upon copy
// seeding is done using previous class.
template<typename Engine, typename Seeder=seeder_t> struct seed_on_copy : Engine{
  typedef Seeder seeder_type;
  seed_on_copy(seeder_type& seeder):Engine(seeder()),s(seeder){}
  seed_on_copy(seed_on_copy const&  o) : Engine(static_cast<Engine const&>(o)), s(o.s)
  { Engine::seed(s()); }
  seeder_type& s;
};

template<std::size_t n_steps=365, std::size_t num_trials=200000, typename Engine=std::mt19937>
struct montecarlo_pricing {
  
  montecarlo_pricing(): years_to_maturity(1.), risk_free_rate(.03), volatility(.2){}
  // we take the engine as parameter to be thread safe : the engine will be thread-private
  template<typename G>
  double step(double s, G& g)const{
    for(int i(0); i != n_steps; ++i)
      { s*= 1. + g()*std::sqrt(years_to_maturity/n_steps) 
	  + risk_free_rate* years_to_maturity/n_steps ; }
    return s;
  }
  // return montecarlo estimate of put/call values
  std::tuple<double, double>
  montecarlo_values(double const strike_price, double const stock_price) const{
    typedef seed_on_copy<Engine> engine_t;
    typename engine_t::seeder_type seeder;
    std::tr1::variate_generator<engine_t, std::tr1::normal_distribution<double> > 
      g(engine_t(seeder), (std::tr1::normal_distribution<double>(0., volatility)));

    std::size_t const nb_threads(std::min(static_cast<std::size_t>(std::max(std::thread::hardware_concurrency(), 1U)), num_trials));
    typedef std::tuple<double, double> result_type; // put and call values, std::pair is so passé :-) 
    std::vector<std::future<result_type> > results;
    std::atomic<std::size_t> remaining_trials(num_trials);
 // there are no std::atomic<> for floating point data, so we cannot have global result (otherwise we could have afforded
 // the small contention at each write
    for(std::size_t i(0); i != nb_threads; ++i){
      results.push_back(std::async([=,&remaining_trials]()-> result_type {
            decltype(g) local_g(g);
            double call_value(0.), put_value(0.);
            while(remaining_trials--){
              double const delta_price(strike_price - this->step(stock_price, local_g));
              *(delta_price >0. ? &put_value : &call_value)+= delta_price;
            }
            return result_type(call_value, put_value);}));
    }
    double call_value(0.), put_value(0.);
    std::for_each(results.begin(), results.end(), [&call_value, &put_value] 
                  (std::future<result_type>& fr){ result_type const& r(fr.get()); call_value-= std::get<0>(r); put_value+= std::get<1>(r); });
    return std::make_tuple(f(call_value), f(put_value));
  }
  // values directly commputed using the Black Scholes formula with constant volatility.
  std::tuple<double, double>
  formula_values(double const strike_price, double const stock_price) const{
    using boost::math::normal;
    normal s;
    double const d1((std::log(stock_price/strike_price)
		     +(risk_free_rate+volatility*volatility/2.)*years_to_maturity)
		    /(volatility*std::sqrt(years_to_maturity)));
    double const d2(d1-volatility*std::sqrt(years_to_maturity));
    double const call(stock_price*cdf(s, d1)-strike_price
		      *std::exp(-risk_free_rate*years_to_maturity)*cdf(s, d2));
    double const put(strike_price*std::exp(-risk_free_rate*years_to_maturity)*cdf(s,-d2)
		     -stock_price*cdf(s, -d1));
    return std::make_tuple(call, put);
  }
  // test main function : computes put & call with monte carlo and forumal.
  void operator()(double const strike_price=110., double const stock_price=100.) const {
    std::tuple<double, double> mc(montecarlo_values(strike_price, stock_price))
      , formula(formula_values(strike_price, stock_price));
    std::cerr<<"delta put:"<<(std::get<1>(mc) - std::get<1>(formula))
	     <<"\tdelta call:"<< (std::get<0>(mc) - std::get<0>(formula))<<std::endl;
    std::cerr<<"mc put:"<<std::get<1>(mc) <<" put:"<<std::get<1>(formula)
	     << "\tmc call:"<< std::get<0>(mc) <<" call:"<<std::get<0>(formula)<<std::endl;
  }
private :
  // just to factor common expression in Monte-Carlo estimate.
  double f(double d) const {return d* std::exp(-risk_free_rate*years_to_maturity) / num_trials ;}
  double const years_to_maturity;
  double const risk_free_rate;
  double const volatility;
};

int main (int argc, char* argv []) {

  montecarlo_pricing<> pricing;
  boost::progress_timer t;
  std::chrono::system_clock::time_point const start (std::chrono::system_clock::now());
  pricing();
  std::chrono::duration<double> const delta(std::chrono::system_clock::now() - start);
  std::cerr<< "wallclock time: "<< delta.count()  <<std::endl;
  return 0;
}
