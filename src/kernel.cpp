
// Author: Manos Papadakis
//[[Rcpp::plugins(cpp11)]]

#define ARMA_64BIT_WORD

#include <RcppArmadillo.h>

#include "Dist.h"
#include "mn.h"

using namespace Rcpp;
using namespace arma;

NumericVector kernel(NumericMatrix X, const double h) {
    const size_t ncl = X.ncol(), nrw = X.nrow();
    mat x(X.begin(), nrw, ncl, false);
    NumericVector Res(ncl);
    colvec res(Res.begin(), ncl, false);
    const double h2 = 2*h*h, k = ( (ncl - 1) * h * sqrt(2 * datum::pi) );
    for (size_t i = 0; i < ncl - 1; ++i) {
        colvec xv(x.begin_col(i), nrw, false);
        long double sv = 0.0;
        for (size_t j = i + 1; j < ncl; ++j) {
            colvec y(x.begin_col(j), nrw, false);
            long double v = exp(-Dist::euclidean<false>(xv, y) / h2);
            sv+=v;
            res[j] += v;
        }
        res[i] += sv;
        res[i] /= k;
    }
    
    res[ncl-1] /= k;
    return Res;
}

RcppExport SEXP Rfast_kernel(SEXP xSEXP, SEXP hSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const double>::type h(hSEXP);
	__result = kernel(x, h);
	return __result;
	END_RCPP
}

#undef ARMA_64BIT_WORD

/////////////////////////////////////////////////////////////////////////////////////
