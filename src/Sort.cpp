// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include "Rfast.h"
#include <vector>
#include <string>
#include "Rfast/parallel.h"

using namespace Rcpp;
using std::remove_if;
using std::string;
using std::vector;
using std::greater;

static SEXP Sort_simple(SEXP x, const bool descend, const bool parallel = false)
{
  SEXP f = PROTECT(Rf_duplicate(x));
  int len = LENGTH(x);
  switch (TYPEOF(x))
  {
  case INTSXP:
  {
    int *F = INTEGER(f);
    descend ? Rfast::sort(F, F + len, std::greater<int>(), parallel) : Rfast::sort(F, F + len, parallel);
    break;
  }
  default:
  {
    double *F = REAL(f);
    descend ? Rfast::sort(F, F + len, std::greater<double>(), parallel) : Rfast::sort(F, F + len, parallel);
    break;
  }
  }
  UNPROTECT(1);
  return f;
}

static NumericVector Sort_na_rm(SEXP x, const bool descend, const bool parallel = false)
{ // na.rm=NA
  NumericVector f(Rf_duplicate(x));
  const int n = remove_if(f.begin(), f.begin() + f.size(), R_IsNA) - f.begin();
  Rfast::sort(f.begin(), f.begin() + n, parallel);
  return f[Range(0, n - 1)];
}

static NumericVector Sort_na_last(SEXP x, const bool descend, const bool parallel = false)
{ // na.rm=T
  NumericVector f(Rf_duplicate(x));
  const int n = remove_if(f.begin(), f.begin() + f.size(), R_IsNA) - f.begin();
  Rfast::sort(f.begin(), f.begin() + n, parallel);
  for (NumericVector::iterator it = f.begin() + n; it != f.end(); ++it)
    *it = NA_REAL;
  return f;
}

SEXP Sort(SEXP x, const bool descend, SEXP na, const bool parallel = false)
{
  if (Rf_isNull(na))
    return Sort_simple(x, descend, parallel);
  else if (R_IsNA(Rf_asReal(na)))
  {
    return Sort_na_rm(x, descend, parallel);
  }
  else
  {
    return Sort_na_last(x, descend, parallel);
  }
  stop("Wrong type of na.last argument.\n");
}

RcppExport SEXP Rfast_Sort(SEXP x, SEXP descendSEXP, SEXP na, SEXP parallelSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  __result = Sort(x, descend, na, parallel);
  return __result;
  END_RCPP
}

//[[Rcpp::export]]
vector<double> Sort_na_first(vector<double> f, const bool descend, const bool parallel = false)
{ // na.rm=F
  const int n = remove_if(f.rbegin(), f.rend(), R_IsNA) - f.rbegin();
  descend ? Rfast::sort(f.end() - n, f.end(), std::greater<double>(), parallel) : Rfast::sort(f.end() - n, f.end(), parallel);
  for (vector<double>::iterator ff = f.begin(); ff != f.begin() + n; ++ff)
    *ff = NA_REAL;
  return f;
}

RcppExport SEXP Rfast_Sort_na_first(SEXP xSEXP, SEXP descendSEXP, SEXP parallelSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<vector<double>>::type x(xSEXP);
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  __result = Sort_na_first(x, descend, parallel);
  return __result;
  END_RCPP
}

//[[Rcpp::export]]
vector<string> sort_string(CharacterVector x, const bool descend, const bool parallel = false)
{
  vector<string> f(x.begin(), x.end());
  descend ? Rfast::sort(f.begin(), f.end(), std::greater<string>(), parallel) : Rfast::sort(f.begin(), f.end(), parallel);
  return f;
}

RcppExport SEXP Rfast_sort_string(SEXP xSEXP, SEXP descendSEXP, SEXP parallelSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<CharacterVector>::type x(xSEXP);
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  __result = sort_string(x, descend, parallel);
  return __result;
  END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////////////////


SEXP stable_sort(SEXP x, const bool descend, const bool parallel = false)
{
  SEXP f = PROTECT(Rf_duplicate(x));
  int len = LENGTH(x);
  switch (TYPEOF(x))
  {
  case INTSXP:
  {
    int *F = INTEGER(f);
    descend ? Rfast::stable_sort(F, F + len, greater<int>(), parallel) : Rfast::stable_sort(F, F + len, parallel);
    break;
  }
  default:
  {
    double *F = REAL(f);
    descend ? Rfast::stable_sort(F, F + len, greater<double>(), parallel) : Rfast::stable_sort(F, F + len, parallel);
    break;
  }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_stable_sort(SEXP x, SEXP descendSEXP, SEXP parallelSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  __result = stable_sort(x, descend, parallel);
  return __result;
  END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////////

static void min_max_neg_pos_helper(vector<int> &x, int &mnn, int &mxn, int &mnp, int &mxp, bool &has_pos, bool &has_neg)
{
  int v;
  for (auto xx = x.begin(); xx != x.end(); ++xx)
  {
    v = *xx;
    if (v < 0)
    {
      has_neg = true;
      if (v < mnn)
        mnn = v;
      else if (v > mxn)
        mxn = v;
    }
    else
    {
      has_pos = true;
      if (v > mxp)
        mxp = v;
      else if (v < mnp)
        mnp = v;
    }
  }
}

//[[Rcpp::export]]
vector<int> sort_int(vector<int> x)
{
  int mnp = INT_MAX, mnn = -1, mxp = 0, mxn = INT_MIN;
  bool has_pos = false, has_neg = false;
  min_max_neg_pos_helper(x, mnn, mxn, mnp, mxp, has_pos, has_neg);
  vector<int> pos, f(x.size()), neg;
  vector<int>::iterator a = x.begin(), F = f.begin();
  if (has_pos)
  {
    pos.resize(mxp - mnp + 1, 0);
  }
  if (has_neg)
  {
    neg.resize(1 - mnn + mxn, 0);
  }
  if (has_pos && has_neg)
  {
    int aa;
    for (; a != x.end(); ++a)
    {
      aa = *a;
      aa < 0 ? ++neg[-aa + mxn] : ++pos[aa - mnp];
    }
  }
  else if (has_pos)
  {
    for (; a != x.end(); ++a)
    {
      ++pos[*a - mnp];
    }
  }
  else
  {
    for (; a != x.end(); ++a)
    {
      ++neg[*a - mxn];
    }
  }

  if (has_neg)
  {
    for (vector<int>::reverse_iterator nr = neg.rbegin(); nr != neg.rend(); ++nr)
    {
      if (*nr != 0)
      {
        int num = -(neg.rend() - nr - 1 - mxn), times = *nr;
        for (int i = 0; i < times; ++i)
        {
          *F++ = num;
        }
      }
    }
  }
  if (has_pos)
  {
    for (a = pos.begin(); a != pos.end(); ++a)
    {
      if (*a != 0)
      {
        int num = a - pos.begin() + mnp, times = *a;
        for (int i = 0; i < times; ++i)
        {
          *F++ = num;
        }
      }
    }
  }
  return f;
}

RcppExport SEXP Rfast_sort_int(SEXP xSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<vector<int>>::type x(xSEXP);
  __result = sort_int(x);
  return __result;
  END_RCPP
}
