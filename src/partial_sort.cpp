//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;

using std::greater;
using std::sort;
using std::nth_element;

IntegerVector partial_sort_index(NumericVector x,const int n,const bool descend){
  IntegerVector ind=seq(1,x.size());
  if(descend){
      auto descend_func = [&](int i,int j){return x[i-1]>x[j-1];};
      nth_element(ind.begin(),ind.begin()+n-1,ind.end(),descend_func);
      sort(ind.begin(),ind.begin()+n,descend_func);
  }else{
      auto descend_func = [&](int i,int j){return x[i-1]<x[j-1];};
      nth_element(ind.begin(),ind.begin()+n-1,ind.end(),descend_func);
      sort(ind.begin(),ind.begin()+n,descend_func);
  }
  return ind;
}

SEXP partial_sort(SEXP x,const int n,const bool descend){
  SEXP f=PROTECT(Rf_duplicate(x));
  int len=LENGTH(x);
  switch(TYPEOF(x)){
	  case INTSXP:{
	    int *F=INTEGER(f);
	    if(descend){
	    	nth_element(F,F+n-1,F+len,greater<int>());
	    	sort(F,F+n,greater<int>());
	    }
	    else {
	    	nth_element(F,F+n-1,F+len);
	    	sort(F,F+n);
	    }
	    break;
	  }
	  default:{
	    double *F=REAL(f);
	    if(descend){
	    	nth_element(F,F+n-1,F+len,greater<double>());
	    	sort(F,F+n,greater<double>());
	    }
	    else {
	    	nth_element(F,F+n-1,F+len);
	    	sort(F,F+n);
	    }
	    break;
	  }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_partial_sort(SEXP x,SEXP nSEXP,SEXP descendSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type n(nSEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    __result = partial_sort(x,n,descend);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_partial_sort_index(SEXP xSEXP,SEXP nSEXP,SEXP descendSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< const int >::type n(nSEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    __result = wrap(partial_sort_index(x,n,descend));
    return __result;
END_RCPP
}