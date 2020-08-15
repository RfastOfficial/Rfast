//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <algorithm>
#include <Rinternals.h>

using namespace Rcpp;
using std::binary_search;
using std::lower_bound;

bool binarysearch(SEXP x,double v){
	bool res;
	switch(TYPEOF(x)){
		case INTSXP:{
			int *start=INTEGER(x);
			res = binary_search(start,start+LENGTH(x),v);
			break;
		}default:{
			double *start=REAL(x);
			res = binary_search(start,start+LENGTH(x),v);
			break;
		}
	}
	return res;
}

int lowerbound(SEXP x,double v){
	int res;
	switch(TYPEOF(x)){
		case INTSXP:{
			int *start=INTEGER(x);
			res = lower_bound(start,start+LENGTH(x),v)-start+1;
			break;
		}default:{
			double *start=REAL(x);
			res = lower_bound(start,start+LENGTH(x),v)-start+1;
			break;
		}
	}
	return res;
}

RcppExport SEXP Rfast_binarysearch(SEXP x,SEXP vSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< double >::type v(vSEXP);
    __result = binarysearch(x,v);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_lowerbound(SEXP x,SEXP vSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< double >::type v(vSEXP);
    __result = lowerbound(x,v);
    return __result;
END_RCPP
}
