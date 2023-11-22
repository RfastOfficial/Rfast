
// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

using std::nth_element;

SEXP col_nth_p(NumericMatrix x, IntegerVector elems, const bool descend, const bool na_rm, const bool index, const unsigned int cores)
{
  const int n = elems.size();
  mat xx(x.begin(), x.nrow(), n, false);
  SEXP F = PROTECT(Rf_allocVector(REALSXP, n));
  double *FF = REAL(F);
  if (index)
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int i = 0; i < n; ++i)
    {
      colvec y = xx.col(i);
      FF[i] = nth_helper_index<colvec>(y, elems[i], descend, na_rm);
    }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int i = 0; i < n; ++i)
    {
      colvec y = xx.col(i);
      FF[i] = nth_helper<colvec>(y, elems[i], descend, na_rm);
    }
  }
  UNPROTECT(1);
  return F;
}

// nth_element
RcppExport SEXP Rfast_col_nth_p(SEXP xSEXP, SEXP ySEXP, SEXP descendSEXP, SEXP na_rmSEXP, SEXP indexSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<IntegerVector>::type y(ySEXP);
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  traits::input_parameter<const bool>::type index(indexSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = col_nth_p(x, y, descend, na_rm, index, cores);
  return __result;
  END_RCPP
}

SEXP row_nth_p(NumericMatrix x, IntegerVector elems, const bool descend, const bool na_rm, const bool index, const unsigned int cores)
{
  const int n = elems.size();
  mat xx(x.begin(), n, x.ncol(), false);
  SEXP F;
  if (index)
  {
    F = PROTECT(Rf_allocVector(INTSXP, n));
    int *FF = INTEGER(F);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int i = 0; i < n; ++i)
    {
      rowvec y = xx.row(i);
      FF[i] = nth_helper_index<rowvec>(y, elems[i], descend, na_rm);
    }
  }
  else
  {
    F = PROTECT(Rf_allocVector(REALSXP, n));
    double *FF = REAL(F);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int i = 0; i < n; ++i)
    {
      rowvec y = xx.row(i);
      FF[i] = nth_helper<rowvec>(y, elems[i], descend, na_rm);
    }
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_row_nth_p(SEXP xSEXP, SEXP ySEXP, SEXP descendSEXP, SEXP na_rmSEXP, SEXP indexSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<IntegerVector>::type y(ySEXP);
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  traits::input_parameter<const bool>::type index(indexSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = row_nth_p(x, y, descend, na_rm, index, cores);
  return __result;
  END_RCPP
}

////////////////////////////////////////////////////////////////////

IntegerMatrix col_order_p(NumericMatrix x, const bool stable, const bool descending, const unsigned int cores)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  IntegerMatrix f(nrw, ncl);
  mat xx(x.begin(), nrw, ncl, false);
  imat ff(f.begin(), nrw, ncl, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
  for (int i = 0; i < ncl; ++i)
  {
    ff.col(i) = Order<icolvec, colvec>(xx.col(i), stable, descending, 1);
  }
  return f;
}

RcppExport SEXP Rfast_col_order_p(SEXP xSEXP, SEXP stableSEXP, SEXP descendingSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<const bool>::type stable(stableSEXP);
  traits::input_parameter<const bool>::type descending(descendingSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = col_order_p(x, stable, descending, cores);
  return __result;
  END_RCPP
}

IntegerMatrix row_order_p(NumericMatrix x, const bool stable, const bool descending, const unsigned int cores)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  IntegerMatrix f(nrw, ncl);
  mat xx(x.begin(), nrw, ncl, false);
  imat ff(f.begin(), nrw, ncl, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
  for (int i = 0; i < nrw; ++i)
  {
    ff.row(i) = Order<irowvec, rowvec>(xx.row(i), stable, descending, 1);
  }
  return f;
}

RcppExport SEXP Rfast_row_order_p(SEXP xSEXP, SEXP stableSEXP, SEXP descendingSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<const bool>::type stable(stableSEXP);
  traits::input_parameter<const bool>::type descending(descendingSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = row_order_p(x, stable, descending, cores);
  return __result;
  END_RCPP
}

////////////////////////////////////////////////////////////////////////////////////

NumericMatrix row_ranks_p(NumericMatrix x, string method, const bool descend, const bool stable, const unsigned int cores)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  NumericMatrix f(nrw, ncl);
  mat xx(x.begin(), nrw, ncl, false);
  mat ff(f.begin(), nrw, ncl, false);
  if (method == "average")
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int i = 0; i < nrw; ++i)
    {
      ff.row(i) = rank_mean<rowvec, rowvec, ivec>(xx.row(i), descend);
    }
  }
  else if (method == "min")
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int i = 0; i < nrw; ++i)
    {
      ff.row(i) = rank_min<rowvec, rowvec, ivec>(xx.row(i), descend);
    }
  }
  else if (method == "max")
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int i = 0; i < nrw; ++i)
    {
      ff.row(i) = rank_max<rowvec, rowvec, ivec>(xx.row(i), descend);
    }
  }
  else if (method == "first")
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int i = 0; i < nrw; ++i)
    {
      ff.row(i) = rank_first<rowvec, rowvec, ivec>(xx.row(i), descend, stable);
    }
  }
  else
    stop("Error. Wrong method.");
  return f;
}

RcppExport SEXP Rfast_row_ranks_p(SEXP xSEXP, SEXP methodSEXP, SEXP descendSEXP, SEXP stableSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<string>::type method(methodSEXP);
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type stable(stableSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = row_ranks_p(x, method, descend, stable, cores);
  return __result;
  END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////

SEXP col_sums_p(SEXP x, const unsigned int cores)
{
  const int n = Rf_ncols(x);
  SEXP F;
  switch (Rfast::Type::type(x))
  {
    case Rfast::Type::Types::REAL:
    {
      F = PROTECT(Rf_allocVector(REALSXP, n));
      double *FF = REAL(F);
      mat xx(REAL(x), Rf_nrows(x), n, false);
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(cores)
      #endif
      for (int i = 0; i < n; i++)
      {
        FF[i] = accu(xx.col(i));
      }
      break;
    }
    default:
    {
      F = PROTECT(Rf_allocVector(INTSXP, n));
      int *FF = INTEGER(F);
      imat xx(INTEGER(x), Rf_nrows(x), n, false);
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(cores)
      #endif
      for (int i = 0; i < n; i++)
      {
        FF[i] = accu(xx.col(i));
      }
      break;
    }
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_col_sums_p(SEXP xSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = col_sums_p(xSEXP, cores);
  return __result;
  END_RCPP
}

SEXP row_sums_p(SEXP x, const unsigned int cores)
{
  const int n = Rf_nrows(x);
  SEXP F;
  switch (Rfast::Type::type(x))
  {
    case Rfast::Type::Types::REAL:
    {
      F = PROTECT(Rf_allocVector(REALSXP, n));
      double *FF = REAL(F);
      mat xx(REAL(x), n, Rf_ncols(x), false);
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(cores)
      #endif
      for (int i = 0; i < n; i++)
      {
        FF[i] = accu(xx.row(i));
      }
      break;
    }
    default:
    {
      F = PROTECT(Rf_allocVector(INTSXP, n));
      int *FF = INTEGER(F);
      imat xx(INTEGER(x), n, Rf_ncols(x), false);
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(cores)
      #endif
      for (int i = 0; i < n; i++)
      {
        FF[i] = accu(xx.row(i));
      }
      break;
    }
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_row_sums_p(SEXP xSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = row_sums_p(xSEXP, cores);
  return __result;
  END_RCPP
}

//////////////////////////////////////////////////////////////

SEXP col_all_p(LogicalMatrix x, const unsigned int cores)
{
  const int n = x.ncol();
  SEXP f = PROTECT(Rf_allocVector(LGLSXP, n));
  imat xx(x.begin(), x.nrow(), n, false);
  int *ff = LOGICAL(x);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
  for (int i = 0; i < n; ++i)
  {
    ff[i] = all(xx.col(i));
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_col_all_p(SEXP xSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<LogicalMatrix>::type x(xSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = col_all_p(x, cores);
  return __result;
  END_RCPP
}

SEXP row_all_p(LogicalMatrix x, const unsigned int cores)
{
  const int n = x.nrow();
  SEXP f = PROTECT(Rf_allocVector(LGLSXP, n));
  imat xx(x.begin(), n, x.ncol(), false);
  int *ff = LOGICAL(x);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
  for (int i = 0; i < n; ++i)
  {
    ff[i] = all(xx.row(i));
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_row_all_p(SEXP xSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<LogicalMatrix>::type x(xSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = row_all_p(x, cores);
  return __result;
  END_RCPP
}

//////////////////////////////////////////////////////////////////////

IntegerVector col_count_values_p(NumericMatrix x, NumericVector values, const unsigned int cores)
{
  const int n = values.size(), p = x.nrow();
  IntegerVector f(n);
  mat xx(x.begin(), p, n, false);
  ivec ff(f.begin(), n, false);
  colvec vv(values.begin(), n, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
  for (int i = 0; i < n; ++i)
  {
    ff[i] = count_value_helper<colvec, double>(xx.col(i), vv[i]);
  }
  return f;
}

RcppExport SEXP Rfast_col_count_values_p(SEXP xSEXP, SEXP valuesSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<NumericVector>::type values(valuesSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = col_count_values_p(x, values, cores);
  return __result;
  END_RCPP
}

IntegerVector row_count_values_p(NumericMatrix x, NumericVector values, const unsigned int cores)
{
  const int n = values.size(), p = x.nrow();
  IntegerVector f(n);
  mat xx(x.begin(), p, n, false);
  ivec ff(f.begin(), n, false);
  colvec vv(values.begin(), n, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
  for (int i = 0; i < n; ++i)
  {
    ff[i] = count_value_helper<rowvec, double>(xx.row(i), vv[i]);
  }
  return f;
}

RcppExport SEXP Rfast_row_count_values_p(SEXP xSEXP, SEXP valuesSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<NumericVector>::type values(valuesSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = row_count_values_p(x, values, cores);
  return __result;
  END_RCPP
}

//////////////////////////////////////////////////////////////////////////