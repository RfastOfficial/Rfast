//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <algorithm>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using std::greater;

//[[Rcpp::export]]
NumericMatrix sort_mat(NumericMatrix x,const bool descend,const bool by_row,const bool stable,const bool parallel, const unsigned int cores){
	return by_row ? Rfast::rowSort(x,descend,stable,parallel,cores) : Rfast::colSort(x,descend,stable,parallel,cores);
}


// sort_mat
RcppExport SEXP Rfast_sort_mat(SEXP xSEXP,SEXP descendSEXP,SEXP by_rowSEXP,SEXP stableSEXP,SEXP parallelSEXP,SEXP coresSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const bool >::type descend(descendSEXP);
    traits::input_parameter< const bool >::type by_row(by_rowSEXP);
    traits::input_parameter< const bool >::type stable(stableSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    traits::input_parameter< const unsigned int >::type cores(coresSEXP);
    if(Rf_isMatrix(xSEXP)){
        __result = sort_mat(NumericMatrix(xSEXP),descend,by_row,stable,parallel,cores);
    }else if(Rf_isNewList(xSEXP)){
        __result = Rfast::colSort(DataFrame(xSEXP),descend,stable,parallel,cores);
    }
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////////////