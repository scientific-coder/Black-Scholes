// Copyright Bernard Hugueney, Licence GPL v3
#include <tr1/random>
#include <numeric>
#include <utility>
#include <boost/progress.hpp>
#include <boost/math/distributions/normal.hpp>
#include <iostream>
#include <time.h>

// g++ openmp-montecarlo.cxx  -o openmp-montecarlo -O4 -march=native -fopenmp -lgomp
// ./openmp-montecarlo


// returns a new seed upon each call
struct seeder_t{
  unsigned int operator()() const{ 
    unsigned long int my_c; 
#pragma omp critical
    { my_c= ++c; }
    return static_cast<unsigned long int>(time(0))+ my_c;
  }
  static unsigned long int c;
};
unsigned long int seeder_t::c=0;
// each thread will get it's own private engine, and each each will have to be re seeded upon copy
// seeding is done using previous class.
template<typename Engine, typename Seeder=seeder_t> struct seed_on_copy : Engine{
  typedef Seeder seeder_type;
  seed_on_copy(seeder_type& seeder):Engine(seeder()),s(seeder){}
  seed_on_copy(seed_on_copy const&  o) : Engine(static_cast<Engine const&>(o)), s(o.s)
  { Engine::seed(s()); }
  seeder_type& s;
};

template<std::size_t n_steps=365, std::size_t num_trials=20000, typename Engine=std::tr1::mt19937>
struct edo_stoch {

  edo_stoch(): years_to_maturity(1.), risk_free_rate(.03), volatility(.2){}
  // we take the engine as parameter to be thread safe : the engine will be thread-private
  template<typename G>
  double step(double s, G& g)const{
    for(int i(1); i != n_steps; ++i)
      { s*= 1. + g()*std::sqrt(years_to_maturity/n_steps) 
	  + risk_free_rate* years_to_maturity/n_steps ; }
    return s;
  }
  // return montecarlo estimate of put/call values
  std::pair<double, double>
  montecarlo_values(double const strike_price, double const stock_price) const{
    double put_value(0.), call_value(0.);
    typedef seed_on_copy<Engine> engine_t;
    typename engine_t::seeder_type seeder;
    std::tr1::variate_generator<engine_t, std::tr1::normal_distribution<double> > 
      g(engine_t(seeder), (std::tr1::normal_distribution<double>(0., volatility)));
#pragma omp parallel for reduction(+: call_value, put_value) firstprivate(g)
      for(std::size_t i=0; i < num_trials; ++i){
	double const delta_price(strike_price - step(stock_price, g));
	*(delta_price >0. ? &put_value : &call_value)+= delta_price;
      }

      return std::make_pair(f(-call_value), f(put_value));
  }
  // values directly commputed using the Black Scholes formula with constant volatility.
  std::pair<double, double>
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
    return std::make_pair(call, put);
  }
  // test main function : computes put & call with monte carlo and forumal.
  void operator()(double const strike_price=110., double const stock_price=100.) const {
    std::pair<double, double> mc(montecarlo_values(strike_price, stock_price))
      , formula(formula_values(strike_price, stock_price));
    std::cerr<<"delta put:"<<(mc.second-formula.second)
	     <<"\tdelta call:"<< (mc.first-formula.first)<<std::endl;
    std::cerr<<"mc put:"<<mc.second<<" put:"<<formula.second
	     << "\tmc call:"<< mc.first<<" call:"<<formula.first<<std::endl;
  }
private :
  // just to factor common expression in Monte-Carlo estimate.
  double f(double d) const {return d* std::exp(-risk_free_rate*years_to_maturity) / num_trials ;}
  double const years_to_maturity;
  double const risk_free_rate;
  double const volatility;
};

int main(int argc, char* argv[]){
  edo_stoch<> s;
  boost::progress_timer t;
  s();
  return 0;
}
