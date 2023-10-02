// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include "Rfast/parallel.h"

using namespace Rcpp;
using namespace std;

using std::greater;
using Rfast::nth_element;
using Rfast::sort;

IntegerVector partial_sort_index(NumericVector x, const int n, const bool descend, const bool parallel = false)
{
  IntegerVector ind = seq(1, x.size());
  if (descend)
  {
    auto descend_func = [&](int i, int j)
    { return x[i - 1] > x[j - 1]; };
    nth_element(ind.begin(), ind.begin() + n - 1, ind.end(), descend_func, parallel);
    sort(ind.begin(), ind.begin() + n, descend_func, parallel);
  }
  else
  {
    auto descend_func = [&](int i, int j)
    { return x[i - 1] < x[j - 1]; };
    nth_element(ind.begin(), ind.begin() + n - 1, ind.end(), descend_func, parallel);
    sort(ind.begin(), ind.begin() + n, descend_func, parallel);
  }
  return ind;
}

SEXP partial_sort(SEXP x, const int n, const bool descend, const bool parallel = false)
{
  SEXP f = PROTECT(Rf_duplicate(x));
  int len = LENGTH(x);
  switch (TYPEOF(x))
  {
  case INTSXP:
  {
    int *F = INTEGER(f);
    if (descend)
    {
      nth_element(F, F + n - 1, F + len, greater<int>(), parallel);
      sort(F, F + n, greater<int>(), parallel);
    }
    else
    {
      nth_element(F, F + n - 1, F + len, parallel);
      sort(F, F + n, parallel);
    }
    break;
  }
  default:
  {
    double *F = REAL(f);
    if (descend)
    {
      nth_element(F, F + n - 1, F + len, greater<double>(), parallel);
      sort(F, F + n, greater<double>(), parallel);
    }
    else
    {
      nth_element(F, F + n - 1, F + len, parallel);
      sort(F, F + n, parallel);
    }
    break;
  }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_partial_sort(SEXP x, SEXP nSEXP, SEXP descendSEXP, SEXP parallelSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const int>::type n(nSEXP);
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  __result = partial_sort(x, n, descend, parallel);
  return __result;
  END_RCPP
}

RcppExport SEXP Rfast_partial_sort_index(SEXP xSEXP, SEXP nSEXP, SEXP descendSEXP, SEXP parallelSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type x(xSEXP);
  traits::input_parameter<const int>::type n(nSEXP);
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  __result = partial_sort_index(x, n, descend, parallel);
  return __result;
  END_RCPP
}