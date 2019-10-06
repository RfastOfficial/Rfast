//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <algorithm>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using std::greater;

//[[Rcpp::export]]
NumericMatrix sort_mat(NumericMatrix x,const bool descend,const bool by_row,const bool stable,const bool parallel){
	return by_row ? Rfast::matrix::rowSort(x,descend,stable,parallel) : Rfast::matrix::colSort(x,descend,stable,parallel);
}

// sort_mat
RcppExport SEXP Rfast_sort_mat(SEXP xSEXP,SEXP descendSEXP,SEXP by_rowSEXP,SEXP stableSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    traits::input_parameter< const bool >::type by_row(by_rowSEXP);
    traits::input_parameter< const bool >::type stable(stableSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = sort_mat(x,descend,by_row,stable,parallel);
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////////////