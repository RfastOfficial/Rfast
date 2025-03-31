//Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <Rinternals.h>
#include <R.h>
#include "mn.h"

using namespace Rcpp;
using std::tgamma;

SEXP Trigamma(SEXP x){
  int n=LENGTH(x);
  SEXP f=PROTECT(Rf_duplicate(x));
  switch(TYPEOF(x)){
    case REALSXP:{
      double *start_f=REAL(f),*start_x=REAL(x),*end_x=start_x+n;
      for(;start_x!=end_x;++start_x,++start_f)
        *start_f=trigamma(*start_x);
      break;
    }
    default:{
      int *start_f=INTEGER(f),*start_x=INTEGER(x),*end_x=start_x+n;
      for(;start_x!=end_x;++start_x,++start_f)
        *start_f=trigamma(*start_x);
      break;
    }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_Trigamma(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = Trigamma(x);
    return __result;
END_RCPP
}


SEXP Digamma(SEXP x){
  int n=LENGTH(x);
  SEXP f=PROTECT(Rf_duplicate(x));
  switch(TYPEOF(x)){
    case REALSXP:{
      double *start_f=REAL(f),*start_x=REAL(x),*end_x=start_x+n;
      for(;start_x!=end_x;++start_x,++start_f)
        *start_f=digamma(*start_x);
      break;
    }
    default:{
      int *start_f=INTEGER(f),*start_x=INTEGER(x),*end_x=start_x+n;
      for(;start_x!=end_x;++start_x,++start_f)
        *start_f=digamma(*start_x);
      break;
    }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_Digamma(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = Digamma(x);
    return __result;
END_RCPP
}

SEXP Lgamma(SEXP x){
  int n=LENGTH(x);
  SEXP f=PROTECT(Rf_duplicate(x));
  switch(TYPEOF(x)){
    case REALSXP:{
      double *start_f=REAL(f),*start_x=REAL(x),*end_x=start_x+n;
      for(;start_x!=end_x;++start_x,++start_f)
        *start_f=lgamma(*start_x);
      break;
    }
    default:{
      int *start_f=INTEGER(f),*start_x=INTEGER(x),*end_x=start_x+n;
      for(;start_x!=end_x;++start_x,++start_f)
        *start_f=lgamma(*start_x);
      break;
    }
  }
  UNPROTECT_PTR(f);
  return f;
}

RcppExport SEXP Rfast_Lgamma(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = Lgamma(x);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////////////////////

SEXP Choose(SEXP n,const int k){
  const int tgk_1=tgamma(k+1);
  int len=LENGTH(n);
  SEXP f=PROTECT(Rf_allocVector(REALSXP,len));
  double *start_f=REAL(f);
  switch(TYPEOF(n)){
    case INTSXP:{
      int *start=INTEGER(n),*end=start+len;
      for(;start!=end;++start,++start_f)
        *start_f=tgamma(*start+1)/(tgk_1*tgamma(*start-k+1));
      break;
    }
    default:{
      double *start=REAL(n),*end=start+len;
      for(;start!=end;++start,++start_f)
        *start_f=tgamma(*start+1)/(tgk_1*tgamma(*start-k+1));
      break;
    }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_Choose(SEXP x,SEXP kSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type k(kSEXP);
    __result = Choose(x,k);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////////////////////


using std::lgamma;

SEXP Lchoose(SEXP n,const int k){
  const double lgk_1=lgamma(k+1);
  int len=LENGTH(n);
  SEXP f=PROTECT(Rf_allocVector(REALSXP,len));
  double *start_f=REAL(f);
  switch(TYPEOF(n)){
  case INTSXP:{
    int *start=INTEGER(n),*end=start+len;
    for(;start!=end;++start,++start_f)
      *start_f=lgamma(*start+1)-lgk_1-lgamma(*start-k+1);
    break;
  }
  default:{
    double *start=REAL(n),*end=start+len;
    for(;start!=end;++start,++start_f)
      *start_f=lgamma(*start+1)-lgk_1-lgamma(*start-k+1);
    break;
  }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_Lchoose(SEXP x,SEXP kSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type k(kSEXP);
    __result = Lchoose(x,k);
    return __result;
END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////

static void combn_mat(arma::vec& vals, const int n, const unsigned int start_idx, 
    std::vector<double>& combn_data, double*& combn_col) {
  if (!n) {
    for (unsigned int i = 0; i < combn_data.size(); ++i) {
      *combn_col++ = combn_data[i];
    }
    return; 
  }
  for (unsigned int i = start_idx; i <= (vals.size() - n); ++i) {
    combn_data.at(combn_data.size() - n) = vals(i);
    combn_mat(vals, n - 1, i + 1, combn_data, combn_col);
  }
}

static void combn_list(arma::vec& vals, const int n, const int start_idx, 
    std::vector<double>& combn_data, int& combn_col,
    Rcpp::List& combn_ds) {
  if (!n) {
    std::vector<double> tmp_combn_data(combn_data.size());
    for (size_t i = 0; i < combn_data.size(); ++i) {
      tmp_combn_data[i] = combn_data[i];
    }
    combn_ds[combn_col++] = tmp_combn_data;
    return;
  }
  for (unsigned int i = start_idx; i <= (vals.size() - n); ++i) {
    combn_data[combn_data.size() - n] = vals[i];
    combn_list(vals, n - 1, i + 1, combn_data, combn_col, combn_ds);
  }
}

SEXP find_combn(arma::vec vals, const int n, const bool ret_mat = true) {
  const unsigned int nrows = n;
  const unsigned int ncols = std::round(R::choose(vals.size(), n));
  std::vector<double> combn_data(nrows);
  const unsigned int start_idx = 0;
  SEXP combn_ds;
  if (ret_mat) {
    static double* combn_col;
    combn_ds = PROTECT(Rf_allocMatrix(REALSXP, nrows, ncols));
    combn_col = REAL(combn_ds); combn_mat(vals, n, start_idx, combn_data, combn_col);
    UNPROTECT(1);
  }
  else {
    static int combn_col;
    Rcpp::List combn_ds_tmp(ncols);
    combn_col = 0; combn_list(vals, n, start_idx, combn_data, combn_col, combn_ds_tmp);
    combn_ds = combn_ds_tmp;
  }
  return combn_ds;
}

RcppExport SEXP Rfast_comb_n(SEXP dataSEXP,SEXP nSEXP,SEXP simplifySEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< arma::vec >::type data(dataSEXP);
    traits::input_parameter< const int >::type n(nSEXP);
    traits::input_parameter< const bool >::type simplify(simplifySEXP);
    __result = find_combn(data,n,simplify);
    return __result;
END_RCPP
}
