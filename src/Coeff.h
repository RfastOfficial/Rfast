

#ifndef COEFF_H
#define COEFF_H

#include <RcppArmadillo.h>
#include <string>

using Rcpp::NumericMatrix;
using std::string;
using namespace arma;
using namespace Rcpp;

namespace Coeff
{
    template<bool Sqrt>
    double bhattacharyya(colvec x, colvec y){
        return Sqrt ? sum(sqrt(x % y)) : dot(x,y);
    }
}

#endif