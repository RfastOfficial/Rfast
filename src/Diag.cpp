
//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <Rinternals.h>
#include <R.h>

using namespace Rcpp;

NumericMatrix diag_matrix_fill_scalar(const int len,const double v){
  SEXP f=PROTECT(Rf_allocMatrix(REALSXP,len,len));
  double *ff=REAL(f),*endf=ff+len*len;
  for(;ff!=endf;++ff){
    *ff=0;
  }
  NumericMatrix x(f);
  x.fill_diag(v);
  UNPROTECT(1);
  return x;
}

SEXP diag_matrix_fill_vec(const int len,SEXP v){
  SEXP f=PROTECT(Rf_allocMatrix(TYPEOF(v),len,len));
  const int len_1=len+1;
  switch(TYPEOF(v)){
  	case REALSXP:{
  	  double *ff=REAL(f),*vv=REAL(v),*endf=ff+len*len;
  	  *ff++=*vv++;
  	  for(int i=1;ff!=endf;++ff,++i){
  	      i==len_1 ? *ff=*vv++,i=0 : *ff=0;
      }
      break;
  	}
  	default:{
  	  int *ff=INTEGER(f),*vv=INTEGER(v),*endf=ff+len*len;
  	  *ff++=*vv++;
  	  for(int i=1;ff!=endf;++ff,++i){
  	      i==len_1 ? *ff=*vv++,i=0 : *ff=0;
      }
  	  break;
    }
  }
  UNPROTECT(1);
  return f;
}

NumericMatrix diag_fill_scalar(NumericMatrix x,const double v){
  NumericMatrix y=clone(x);
  y.fill_diag(v);
  return y;
}

SEXP diag_fill_vec(SEXP x,SEXP v){
  SEXP f=PROTECT(Rf_duplicate(x));
  const int len=Rf_ncols(x),len_1=len+1;
  switch(TYPEOF(x)){
  	case REALSXP:{
  	  double *ff=REAL(f),*vv=REAL(v),*endv=vv+LENGTH(v);
  	  for(;vv!=endv;ff+=len_1,++vv){
  	    *ff=*vv;
  	  }
      break;
  	}
  	default:{
  	  int *ff=INTEGER(f),*vv=INTEGER(v),*endv=vv+LENGTH(v);
  	  for(;vv!=endv;ff+=len_1,++vv){
  	    *ff=*vv;
  	  }
      break;
  	}
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_diag_matrix_fill_scalar(SEXP lenSEXP,SEXP vSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type len(lenSEXP);
    traits::input_parameter< const double >::type v(vSEXP);
    __result = wrap(diag_matrix_fill_scalar(len,v));
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_diag_matrix_fill_vec(SEXP lenSEXP,SEXP v) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type len(lenSEXP);
    __result = diag_matrix_fill_vec(len,v);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_diag_fill_scalar(SEXP xSEXP,SEXP vSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const double >::type v(vSEXP);
    __result = wrap(diag_fill_scalar(x,v));
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_diag_fill_vec(SEXP x,SEXP v) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = diag_fill_vec(x,v);
    return __result;
END_RCPP
}