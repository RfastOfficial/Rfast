

#ifndef COEFF_H
#define COEFF_H

#include <Rcpp.h>
#include <string>

using Rcpp::NumericMatrix;
using std::string;
using namespace arma;

namespace Coeff
{
    inline double bhattacharyya(colvec x, colvec y){
        return sum(sqrt(x % y));
    }
}

#endif