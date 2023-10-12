// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <algorithm>
#include "mn.h"
#include "Rfast/parallel.h"

using namespace Rcpp;
using std::remove_if;

int nth_int(vector<int> x, int elem)
{
  int aa, mx, mn;
  bool has_pos = false, has_neg = false;
  max_neg_pos(&x[0], &x[x.size() - 1] + 1, mx, mn, has_pos, has_neg);
  vector<int> pos, f(x.size()), neg;
  vector<int>::iterator a = x.begin();
  if (has_pos)
  {
    pos.resize(mx + 1, 0);
  }
  if (has_neg)
  {
    neg.resize(1 - mn, 0);
  }
  if (has_pos && has_neg)
  {
    for (; a != x.end(); ++a)
    {
      aa = *a;
      aa < 0 ? ++neg[-aa] : ++pos[aa];
    }
  }
  else if (has_pos)
  {
    for (; a != x.end(); ++a)
    {
      aa = *a;
      ++pos[aa];
    }
  }
  else
  {
    for (; a != x.end(); ++a)
    {
      aa = *a;
      ++neg[-aa];
    }
  }

  --elem;
  int res = 0, num = 0;
  if (has_neg)
  {
    for (vector<int>::reverse_iterator nr = neg.rbegin(); nr != neg.rend(); ++nr)
    {
      if (*nr != 0)
      {
        num += *nr;
        if (elem > num)
          res = nr - neg.rbegin() + 1;
      }
    }
  }
  if (has_pos)
  {
    for (a = pos.begin(); a != pos.end(); ++a)
    {
      if (*a != 0)
      {
        num += *a;
        if (elem > num)
          res = a - pos.begin() + 1;
      }
    }
  }
  return res;
}

// nth
RcppExport SEXP Rfast_nth(SEXP xSEXP, SEXP elemSEXP, SEXP num_of_nthsSEXP, SEXP descendSEXP, SEXP na_rmSEXP, SEXP indexSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  // traits::input_parameter< NumericVector >::type X(xSEXP);
  traits::input_parameter<const int>::type elem(elemSEXP);
  traits::input_parameter<const int>::type num_of_nths(num_of_nthsSEXP);
  traits::input_parameter<const bool>::type descend(descendSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  traits::input_parameter<const bool>::type index(indexSEXP);

  NumericVector x = clone(xSEXP);

  if (num_of_nths > 1)
  {
    colvec y(x.begin(), x.size(), false);
    if (index)
    {
      __result = nth_helper_index_n_elems<colvec>(y, elem, descend, na_rm);
    }
    else
    {
      __result = nth_helper_n_elems<colvec>(y, elem, descend, na_rm);
    }
  }
  else
  {
    if (index)
    {
      __result = nth_helper_index<NumericVector>(x, elem, descend, na_rm);
    }
    else
    {
      __result = nth_helper<NumericVector>(x, elem, descend, na_rm);
    }
  }
  return __result;
  END_RCPP
}

// nth_int
RcppExport SEXP Rfast_nth_int(SEXP xSEXP, SEXP elemSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<vector<int>>::type x(xSEXP);
  traits::input_parameter<int>::type elem(elemSEXP);
  __result = nth_int(x, elem);
  return __result;
  END_RCPP
}

//////////////////////////////////////////////////////////////////////////////

double med(SEXP x, const bool na_rm)
{
  double s = 0;
  switch (TYPEOF(x))
  {
  case REALSXP:
  {
    NumericVector xx(Rf_duplicate(x));
    s = na_rm ? med_helper<NumericVector>(xx.begin(), xx.begin() + (std::remove_if(xx.begin(), xx.end(), R_IsNA) - xx.begin())) : med_helper<NumericVector>(xx.begin(), xx.end());
    break;
  }
  case INTSXP:
  {
    IntegerVector xx(Rf_duplicate(x));
    s = na_rm ? med_helper<IntegerVector>(xx.begin(), xx.begin() + (std::remove_if(xx.begin(), xx.end(), R_IsNA) - xx.begin())) : med_helper<IntegerVector>(xx.begin(), xx.end());
    break;
  }
  default:
    stop("Error: Unknown type.\n");
    break;
  }
  return s;
}

// returns median of a vector
RcppExport SEXP Rfast_med(SEXP x, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = med(x, na_rm);
  return __result;
  END_RCPP
}

///////////////////////////////////////////////////////////////////////////////

SEXP min_max(SEXP x, bool index = false)
{
  SEXP F;
  double *xx = REAL(x), *end = xx + LENGTH(x);
  double xxx;
  if (index)
  {
    F = PROTECT(Rf_allocVector(INTSXP, 2));
    int *f = INTEGER(F), min_i = 0, max_i = 0;
    double *bg = xx;
    for (xx++; xx != end; ++xx)
    {
      xxx = *xx;
      if (xxx > bg[max_i])
      {
        max_i = xx - bg;
      }
      else if (xxx < bg[min_i])
      {
        min_i = xx - bg;
      }
    }
    *f = min_i + 1;
    f[1] = max_i + 1;
    UNPROTECT(1);
    return F;
  }
  F = PROTECT(Rf_allocVector(REALSXP, 2));
  double *f = REAL(F), min = *xx, max = min;
  for (xx++; xx != end; ++xx)
  {
    xxx = *xx;
    if (xxx > max)
      max = xxx;
    else if (xxx < min)
      min = xxx;
  }
  *f = min;
  f[1] = max;
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_min_max(SEXP x, SEXP indexSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<bool>::type index(indexSEXP);
  __result = min_max(x, index);
  return __result;
  END_RCPP
}

SEXP min_max_perc(SEXP x)
{
  const int n = LENGTH(x);
  SEXP f = Rf_allocVector(REALSXP, 4);
  double *start = REAL(x), *end = start + n, mx, mn, pos = 0, xx, *FF = REAL(f);
  mn = mx = *start;
  for (; start != end; ++start)
  {
    xx = *start;
    if (xx > 0)
      pos++;
    if (mn > xx)
      mn = xx;
    else if (mx < xx)
      mx = xx;
  }
  *FF = mn;
  FF[1] = mx;
  FF[3] = (pos / n) * 100.0;
  FF[2] = 100.0 - FF[3];
  return f;
}

RcppExport SEXP Rfast_min_max_perc(SEXP x)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  __result = min_max_perc(x);
  return __result;
  END_RCPP
}

////////////////////////////////////////////////////////////////////

SEXP pmax_simple(SEXP x, SEXP y)
{
  SEXP f = (Rf_isMatrix(x) && Rf_isMatrix(y)) ? PROTECT(Rf_allocMatrix(REALSXP, Rf_nrows(x), Rf_ncols(x))) : PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
  double *startx = REAL(x), *end = startx + LENGTH(x), *starty = REAL(y), *startf = REAL(f);
  for (; startx != end; ++startx, ++starty, ++startf)
    *startf = std::max(*startx, *starty);
  UNPROTECT(1);
  return f;
}

SEXP pmax_na_rm(SEXP x, SEXP y)
{
  SEXP f = (Rf_isMatrix(x) && Rf_isMatrix(y)) ? PROTECT(Rf_allocMatrix(REALSXP, Rf_nrows(x), Rf_ncols(x))) : PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
  double *startx = REAL(x), *end = startx + LENGTH(x), *starty = REAL(y), *startf = REAL(f);
  for (; startx != end; ++startx, ++starty, ++startf)
  {
    if (!(R_IsNA(*startx) || R_IsNA(*starty)))
      *startf = std::max(*startx, *starty);
  }
  UNPROTECT(1);
  return f;
}

SEXP pmax(SEXP x, SEXP y, const bool na_rm)
{
  return na_rm ? pmax_na_rm(x, y) : pmax_simple(x, y);
}

RcppExport SEXP Rfast_pmax(SEXP x, SEXP y, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = pmax(x, y, na_rm);
  return __result;
  END_RCPP
}

SEXP pmin_simple(SEXP x, SEXP y)
{
  SEXP f = (Rf_isMatrix(x) && Rf_isMatrix(y)) ? PROTECT(Rf_allocMatrix(REALSXP, Rf_nrows(x), Rf_ncols(x))) : PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
  double *startx = REAL(x), *end = startx + LENGTH(x), *starty = REAL(y), *startf = REAL(f);
  for (; startx != end; ++startx, ++starty, ++startf)
    *startf = std::min(*startx, *starty);
  UNPROTECT(1);
  return f;
}

SEXP pmin_na_rm(SEXP x, SEXP y)
{
  SEXP f = (Rf_isMatrix(x) && Rf_isMatrix(y)) ? PROTECT(Rf_allocMatrix(REALSXP, Rf_nrows(x), Rf_ncols(x))) : PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
  double *startx = REAL(x), *end = startx + LENGTH(x), *starty = REAL(y), *startf = REAL(f);
  for (; startx != end; ++startx, ++starty, ++startf)
  {
    if (!(R_IsNA(*startx) || R_IsNA(*starty)))
      *startf = std::min(*startx, *starty);
  }
  UNPROTECT(1);
  return f;
}

SEXP pmin(SEXP x, SEXP y, const bool na_rm)
{
  return na_rm ? pmin_na_rm(x, y) : pmin_simple(x, y);
}

RcppExport SEXP Rfast_pmin(SEXP x, SEXP y, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = pmin(x, y, na_rm);
  return __result;
  END_RCPP
}

//[[Rcpp::export]]
SEXP pmin_pmax_simple(SEXP x, SEXP y)
{
  SEXP f = PROTECT(Rf_allocMatrix(REALSXP, 2, LENGTH(x)));
  double *startx = REAL(x), *end = startx + LENGTH(x), *starty = REAL(y), *startf = REAL(f);
  for (; startx != end; ++startx, ++starty, startf += 2)
    if (*startx < *starty)
    {
      *startf = *startx;
      startf[1] = *starty;
    }
    else
    {
      *startf = *starty;
      startf[1] = *startx;
    }
  UNPROTECT(1);
  return f;
}

//[[Rcpp::export]]
SEXP pmin_pmax_na_rm(SEXP x, SEXP y)
{
  SEXP f = PROTECT(Rf_allocMatrix(REALSXP, 2, LENGTH(x)));
  double *startx = REAL(x), *end = startx + LENGTH(x), *starty = REAL(y), *startf = REAL(f), vx, vy;
  for (; startx != end; ++startx, ++starty, startf += 2)
  {
    vx = *startx;
    vy = *starty;
    if (!(R_IsNA(vx) || (R_IsNA(vy))))
    {
      if (vx < vy)
      {
        *startf = vx;
        startf[1] = vy;
      }
      else
      {
        *startf = vy;
        startf[1] = vx;
      }
    }
  }
  UNPROTECT(1);
  return f;
}

SEXP pmin_pmax(SEXP x, SEXP y, const bool na_rm)
{
  return na_rm ? pmin_pmax_na_rm(x, y) : pmin_pmax_simple(x, y);
}

RcppExport SEXP Rfast_pmin_pmax(SEXP x, SEXP y, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = pmin_pmax(x, y, na_rm);
  return __result;
  END_RCPP
}

////////////////////////////////////////////////////////////////////

NumericVector min_freq_d(NumericVector x, const int na_rm)
{
  NumericVector xx = clone(x);
  const int n_1 = na_rm ? x.size() : remove_if(xx.begin(), xx.end(), R_IsNA) - xx.begin();
  int i, j = 0;
  std::sort(xx.begin(), xx.begin() + n_1);
  if (!na_rm)
  {
    xx.push_back(0.0);
  }
  int times = 0;
  double v = xx[j], mn_val = 0.0;
  int mn_fr = INT_MAX;
  for (i = 1; i < n_1; ++i)
  {
    if (v != xx[i])
    {
      times = i - j;
      if (times < mn_fr)
      {
        mn_fr = times;
        mn_val = v;
        if (times == 1)
        {
          break;
        }
      }
      j = i;
      v = xx[j];
    }
  }
  return NumericVector::create(_["value"] = mn_val, _["freq"] = mn_fr);
}

RcppExport SEXP Rfast_min_freq_d(SEXP xSEXP, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type x(xSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = min_freq_d(x, na_rm);
  return __result;
  END_RCPP
}

IntegerVector min_freq_i(IntegerVector X, const bool na_rm)
{
  int sz;
  IntegerVector x;
  if (na_rm)
  {
    x = clone(X);
    sz = remove_if(x.begin(), x.end(), R_IsNA) - x.begin();
  }
  else
  {
    x = X;
    sz = x.size();
  }
  int aa, szp = sz, szn = sz, count_neg = 0, count_pos = 0;
  vector<int> f(sz), ff, neg(sz);
  IntegerVector::iterator a = x.begin();
  vector<int>::iterator F = f.begin(), nn = neg.begin(), index;
  for (; a != x.end(); ++a)
  {
    aa = *a;
    if (aa < 0)
    {
      if (-aa >= szn)
      {
        neg.resize(-aa + 1);
        nn = neg.begin();
        szn = neg.size();
      }
      count_neg++;
      *(nn - aa) += 1;
    }
    else
    {
      if (aa >= szp)
      {
        f.resize(aa + 1);
        F = f.begin();
        szp = f.size();
      }
      count_pos++;
      *(F + aa) += 1;
    }
  }
  int val, freq;
  if (!count_neg)
  {
    index = min_element(f.begin(), f.end());
    val = index - f.begin();
    freq = *index;
  }
  else if (!count_pos)
  {
    index = min_element(neg.begin(), neg.end());
    val = index - f.begin();
    freq = *index;
  }
  else
  {
    vector<int>::iterator index_neg = min_element(neg.begin(), neg.end());
    index = min_element(f.begin(), f.end());
    if (*index > *index_neg)
    {
      freq = *index;
      val = index - f.begin();
    }
    else
    {
      freq = *index_neg;
      val = index_neg - neg.begin();
    }
  }
  return IntegerVector::create(_["value"] = val, _["frequency"] = freq);
}

RcppExport SEXP Rfast_min_freq_i(SEXP xSEXP, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<IntegerVector>::type x(xSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = min_freq_i(x, na_rm);
  return __result;
  END_RCPP
}

NumericVector max_freq_d(NumericVector x, const int na_rm)
{
  NumericVector xx = clone(x);
  const int n_1 = na_rm ? x.size() : remove_if(xx.begin(), xx.end(), R_IsNA) - xx.begin();
  int i, j = 0;
  std::sort(xx.begin(), xx.begin() + n_1);
  if (!na_rm)
  {
    xx.push_back(0.0);
  }
  int times = 0;
  double v = xx[j], mx_val = 0.0;
  int mx_fr = 0;
  for (i = 1; i < n_1; ++i)
  {
    if (v != xx[i])
    {
      times = i - j;
      if (times > mx_fr)
      {
        mx_fr = times;
        mx_val = v;
      }
      j = i;
      v = xx[j];
    }
  }
  return NumericVector::create(_["value"] = mx_val, _["freq"] = mx_fr);
}

RcppExport SEXP Rfast_max_freq_d(SEXP xSEXP, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type x(xSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = max_freq_d(x, na_rm);
  return __result;
  END_RCPP
}

IntegerVector max_freq_i(IntegerVector X, const bool na_rm)
{
  int sz;
  IntegerVector x;
  if (na_rm)
  {
    x = clone(X);
    sz = remove_if(x.begin(), x.end(), R_IsNA) - x.begin();
  }
  else
  {
    x = X;
    sz = x.size();
  }
  int aa, szp = sz, szn = sz, count_neg = 0, count_pos = 0;
  vector<int> f(sz), ff, neg(sz);
  IntegerVector::iterator a = x.begin();
  vector<int>::iterator F = f.begin(), nn = neg.begin(), index;
  for (; a != x.end(); ++a)
  {
    aa = *a;
    if (aa < 0)
    {
      if (-aa >= szn)
      {
        neg.resize(-aa + 1);
        nn = neg.begin();
        szn = neg.size();
      }
      count_neg++;
      *(nn - aa) += 1;
    }
    else
    {
      if (aa >= szp)
      {
        f.resize(aa + 1);
        F = f.begin();
        szp = f.size();
      }
      count_pos++;
      *(F + aa) += 1;
    }
  }
  int val, freq;
  if (!count_neg)
  {
    index = max_element(f.begin(), f.end());
    val = index - f.begin();
    freq = *index;
  }
  else if (!count_pos)
  {
    index = max_element(neg.begin(), neg.end());
    val = index - f.begin();
    freq = *index;
  }
  else
  {
    vector<int>::iterator index_neg = max_element(neg.begin(), neg.end());
    index = max_element(f.begin(), f.end());
    if (*index > *index_neg)
    {
      freq = *index;
      val = index - f.begin();
    }
    else
    {
      freq = *index_neg;
      val = index_neg - neg.begin();
    }
  }
  return IntegerVector::create(_["value"] = val, _["frequency"] = freq);
}

RcppExport SEXP Rfast_max_freq_i(SEXP xSEXP, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<IntegerVector>::type x(xSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = max_freq_i(x, na_rm);
  return __result;
  END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////

SEXP Outer(SEXP x, SEXP y, const char oper = '*')
{
  int lenx = LENGTH(x), leny = LENGTH(y);
  SEXP F = PROTECT(Rf_allocMatrix(REALSXP, leny, lenx));
  double *xx = REAL(x), *end = xx + lenx, *f = REAL(F), *yy = REAL(y);

  switch (oper)
  {
  case '*':
    for (; xx != end; ++xx, f += leny)
      myoperator<double, mmult<double>>(f, *xx, yy, leny);
    break;
  case '-':
    for (; xx != end; ++xx, f += leny)
      myoperator<double, mdiff<double>>(f, *xx, yy, leny);
    break;
  case '+':
    for (; xx != end; ++xx, f += leny)
      myoperator<double, madd<double>>(f, *xx, yy, leny);
    break;
  case '/':
    for (; xx != end; ++xx, f += leny)
      myoperator<double, mdiv<double>>(f, *xx, yy, leny);
    break;
  case '^':
    for (; xx != end; ++xx, f += leny)
      myoperator<double, std::pow>(f, *xx, yy, leny);
    break;
  case '%':
    for (; xx != end; ++xx, f += leny)
      myoperator<double, std::fmod>(f, *xx, yy, leny);
    break;
  default:
    stop("Wrong operator.\n");
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_Outer(SEXP x, SEXP y, SEXP operSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const char>::type oper(operSEXP);
  __result = Outer(x, y, oper);
  return __result;
  END_RCPP
}

///////////////////////////////////////////////////////////////////////////////

SEXP Round_simple(SEXP x, int dg)
{
  const int n = LENGTH(x);
  SEXP f = PROTECT(Rf_duplicate(x));
  double *start = REAL(x), *end = start + n, *ff = REAL(f);
  for (; start != end; ++start, ++ff)
    *ff = my_round_gen_simple(*start, dg);
  UNPROTECT(1);
  return f;
}

SEXP Round_na_rm(SEXP x, const int dg)
{
  const int n = LENGTH(x);
  SEXP f = PROTECT(Rf_duplicate(x));
  double *start = REAL(x), *end = start + n, *ff = REAL(f);
  for (; start != end; ++start, ++ff)
    *ff = my_round_gen_na_rm(*start, dg);
  UNPROTECT(1);
  return f;
}

//[[Rcpp::export]]
SEXP Round(SEXP x, const int dg, const bool na_rm)
{
  return na_rm ? Round_simple(x, dg > 15 ? 15 : dg) : Round_na_rm(x, dg > 15 ? 15 : dg);
}

RcppExport SEXP Rfast_Round(SEXP x, SEXP dgSEXP, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const int>::type dg(dgSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = Round(x, dg, na_rm);
  return __result;
  END_RCPP
}

///////////////////////////////////////////////////////////////////////////////
//[[Rcpp::export]]
NumericMatrix squareform_c(NumericVector x)
{
  const int d = my_round(0.5 + sqrt(1 + 8 * x.size()) / 2.0);
  int i, j, s = 0;
  double a;
  NumericMatrix f(d, d);
  for (i = 0; i < d; ++i)
  {
    for (j = i + 1; j < d; ++j, ++s)
    {
      a = x[s];
      f(j, i) = a;
      f(i, j) = a;
    }
  }
  return f;
}

RcppExport SEXP Rfast_squareform_c(SEXP xSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type x(xSEXP);
  __result = squareform_c(x);
  return __result;
  END_RCPP
}

//////////////////////////////////////////////////////////////////////////////

RcppExport SEXP Rfast_symmetric(SEXP xSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  __result = Rfast::is_symmetric(x);
  return __result;
  END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP Rfast_var(SEXP xSEXP, SEXP stdSEXP, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type x(xSEXP);
  traits::input_parameter<const bool>::type std(stdSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  __result = Rfast::var<NumericVector>(x, std, na_rm);
  return __result;
  END_RCPP
}

///////////////////////////////////////////////////////////////////////////////

int count_value(SEXP x, SEXP value)
{
  int s = 0;
  switch (TYPEOF(value))
  {
  case REALSXP:
    s = count_value_helper<NumericVector, double>(NumericVector(x), Rf_asReal(value));
    break;
  case INTSXP:
    s = count_value_helper<IntegerVector, int>(IntegerVector(x), Rf_asInteger(value));
    break;
  case STRSXP:
    s = count_value_helper<vector<string>, string>(as<vector<string>>(x), as<string>(value));
    break;
  default:
    stop("Error: Unknown type of argument value.\n");
    break;
  }
  return s;
}

RcppExport SEXP Rfast_count_value(SEXP x, SEXP value)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  __result = count_value(x, value);
  return __result;
  END_RCPP
}

////////////////////////////////////////////////////////////////////////////////

SEXP Log(SEXP x)
{
  const int nrow = Rf_nrows(x), ncol = Rf_ncols(x), n = nrow * ncol;
  SEXP f;
  switch (TYPEOF(x))
  {
  case REALSXP:
  {
    f = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
    double *start_f = REAL(f), *start_x = REAL(x), *end_x = start_x + n;
    for (; start_x != end_x; ++start_x, ++start_f)
      *start_f = std::log(*start_x);
    break;
  }
  default:
  {
    f = PROTECT(Rf_allocMatrix(INTSXP, nrow, ncol));
    int *start_f = INTEGER(f), *start_x = INTEGER(x), *end_x = start_x + n;
    for (; start_x != end_x; ++start_x, ++start_f)
      *start_f = std::log(*start_x);
    break;
  }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_Log(SEXP x)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  __result = Log(x);
  return __result;
  END_RCPP
}

//////////////////////////////////////////////////////////////////////////

using std::lgamma;

SEXP Lbeta(SEXP x, SEXP y)
{
  int n = LENGTH(x);
  SEXP f = PROTECT(Rf_duplicate(x));
  switch (TYPEOF(x))
  {
  case REALSXP:
  {
    double *start_f = REAL(f), *start_x = REAL(x), *start_y = REAL(y), *end_x = start_x + n, X, Y;
    for (; start_x != end_x; ++start_x, ++start_y, ++start_f)
    {
      X = *start_x;
      Y = *start_y;
      *start_f = lgamma(X) + lgamma(Y) - lgamma(X + Y);
    }
    break;
  }
  default:
  {
    int *start_f = INTEGER(f), *start_x = INTEGER(x), *start_y = INTEGER(y), *end_x = start_x + n, X, Y;
    for (; start_x != end_x; ++start_x, ++start_y, ++start_f)
    {
      X = *start_x;
      Y = *start_y;
      *start_f = lgamma(X) + lgamma(Y) - lgamma(X + Y);
    }
    break;
  }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_Lbeta(SEXP x, SEXP y)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  __result = Lbeta(x, y);
  return __result;
  END_RCPP
}

//////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
IntegerVector Match(NumericVector x, NumericVector key)
{
  return match(x, key);
}

RcppExport SEXP Rfast_Match(SEXP xSEXP, SEXP keySEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type x(xSEXP);
  traits::input_parameter<NumericVector>::type key(keySEXP);
  __result = Match(x, key);
  return __result;
  END_RCPP
}

////////////////////////////////////////////////////////////////////////////////////

typedef Rboolean (*R_Function)(SEXP);

template <R_Function func>
void which_is_helper(DataFrame &x, vector<int> &P)
{
  DataFrame::iterator xx = x.begin();
  for (int i = 1; xx != x.end(); ++xx, ++i)
  {
    if (func(*xx))
    {
      P.push_back(i);
    }
  }
}

//[[Rcpp::export]]
vector<int> which_is(DataFrame x, const string method)
{
  vector<int> P;
  if (method == "logical")
  {
    which_is_helper<Rf_isLogical>(x, P);
  }
  else if (method == "integer")
  {
    which_is_helper<Rf_isInteger>(x, P);
  }
  else if (method == "factor")
  {
    which_is_helper<Rf_isFactor>(x, P);
  }
  else if (method == "numeric")
  {
    which_is_helper<Rf_isNumeric>(x, P);
  }
  return P;
}

RcppExport SEXP Rfast_which_is(SEXP xSEXP, SEXP methodSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<DataFrame>::type x(xSEXP);
  traits::input_parameter<const string>::type method(methodSEXP);
  __result = which_is(x, method);
  return __result;
  END_RCPP
}

///////////////////////////////////////////////////////////////////////////

RcppExport SEXP Rfast_mad2(SEXP xSEXP, SEXP methodSEXP, SEXP na_rmSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  string method = as<string>(methodSEXP);
  traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
  // if method is median then copy the vector because median changes the memory
  traits::input_parameter<NumericVector>::type x(method == "median" ? Rf_duplicate(xSEXP) : xSEXP);
  __result = Rfast::mad<NumericVector>(x, method, na_rm);
  return __result;
  END_RCPP
}

///////////////////////////////////////////////////////////////////////////////
