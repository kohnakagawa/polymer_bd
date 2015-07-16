#pragma once
#include <cmath>
#include <boost/random.hpp>

class PRNG{
  const size_t n_rng;
  boost::mt19937* gen;
  boost::uniform_real<> uni_dist; //range is [0,1].
  boost::normal_distribution<> *normal_dist; //mean=0 sd=1.0
  
public:
  PRNG(const size_t seed, const size_t th_numb) : n_rng(th_numb) {
    gen = new boost::mt19937 [th_numb];
    for(size_t i = 0; i < th_numb; i++) gen[i].seed(seed + i * 10);
    normal_dist = new boost::normal_distribution<> (0.0, 1.0);
  }
  PRNG(const size_t th_numb) : n_rng(th_numb) {
    gen = new boost::mt19937 [th_numb];
    const size_t base_seed = (size_t) time(NULL);
    for(size_t i = 0; i < th_numb; i++) gen[i].seed(base_seed + i * 10);
    normal_dist = new boost::normal_distribution<> (0.0, 1.0);
  }
  
  ~PRNG(){
    delete [] gen;
    delete normal_dist;
  }
  
  inline double Uniform(const int tid){
    return uni_dist(gen[tid]);
  }
  inline double Uniform(const int tid, const double up, const double dw){
    return dw + (up - dw) * uni_dist(gen[tid]);
  }

  inline double Normal(const int tid){
    return (*normal_dist)(gen[tid]);
  }
  inline double Normal(const int tid, const double mean, const double sd){
    return mean + (*normal_dist)(gen[tid]) * sqrt(sd);
  }
};
