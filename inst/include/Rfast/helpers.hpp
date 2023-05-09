
#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;
using namespace Rcpp;

inline long long int get_current_nanoseconds(){
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

#endif
