
// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <string>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;

using std::string;

static int proper_size(int nrw, int ncl)
{
  return ncl * (ncl - 1) * 0.5;
}

NumericVector euclidean_vec(NumericMatrix x, const bool sqr)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;
  if (sqr)
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j, ++k)
      {
        f[k] = sum(square(xx.col(j) - xv));
      }
    }
  else
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j, ++k)
      {
        f[k] = std::sqrt(sum(square(xv - xx.col(j))));
      }
    }
  return f;
}

NumericVector manhattan_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum(abs(xv - xx.col(j)));
    }
  }
  return f;
}

NumericVector hellinger_vec(NumericMatrix x, const bool sqr)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  const double p = 1.0 / std::sqrt(2.0);
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;
  if (sqr)
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j, ++k)
      {
        f[k] = sum(square(xv - xx.col(j))) * 0.5;
      }
    }
  else
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j, ++k)
      {
        f[k] = p * std::sqrt(sum(square(xv - xx.col(j))));
      }
    }
  return f;
}

NumericVector max_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw), tmp(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      tmp = abs(xv - xx.col(j));
      f[k] = tmp.at(tmp.index_max());
    }
  }
  return f;
}

NumericVector min_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw), tmp(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      tmp = abs(xx.col(j) - xv);
      f[k] = tmp(tmp.index_min());
    }
  }
  return f;
}

NumericVector minkowski_vec(NumericMatrix x, const double p)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  const double p_1 = 1.0 / p;
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = pow(sum_with<std::pow, colvec>(abs(xv - xx.col(j)), p), p_1);
    }
  }
  return f;
}

NumericVector canberra_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw), absx(nrw);
  mat x_abs = abs(x);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    absx = x_abs.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum(abs(xv - xx.col(j)) / (absx + x_abs.col(j)));
    }
  }
  return f;
}

NumericVector total_variation_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = 0.5 * sum(abs(xv - xx.col(j)));
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericVector kullback_leibler_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  NumericVector f(proper_size(nrw, ncl));
  mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
  colvec xv(nrw), log_xv(nrw);
  int i, j, k = 0;
  fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    log_xv = log_xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum((xv - xx.col(j)) % (log_xv - log_xx.col(j)));
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericVector bhattacharyya_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum(sqrt(xv % xx.col(j)));
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericVector itakura_saito_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  NumericVector f(proper_size(nrw, ncl));
  mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
  colvec xv(nrw), log_xv(nrw);
  int i, j, k = 0;
  fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    log_xv = log_xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum(xv / xx.col(j) - (log_xv - log_xx.col(j)) - 1);
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericVector jensen_shannon_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  NumericVector f(proper_size(nrw, ncl));
  mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
  colvec xv(nrw), log_xv(nrw);
  constexpr double log2 = std::log(2);
  int i, j, k = 0;
  fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    log_xv = log_xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum_with_condition<double, check_if_is_finite, colvec>((xv + xx.col(j)) % (log2 - arma::log(xv + xx.col(j))) + xv % log_xv + xx.col(j) % log_xx.col(j));
    }
  }
  return f;
}

NumericVector cosine_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  NumericVector f(proper_size(nrw, ncl));
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw), norm_x = euclidean_norm(xx);
  int i, j, k = 0;

  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    double normx = norm_x[i];
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = dot(xv, xx.col(j)) / (normx * norm_x[j]);
    }
  }
  return f;
}

NumericVector soergel_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum(abs(xv - xx.col(j))) / sum_max_elems(xv, xx.col(j));
    }
  }
  return f;
}

NumericVector chi_square_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum(square(xv - xx.col(j)) / (xv + xx.col(j)));
    }
  }
  return f;
}

NumericVector sorensen_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j, k = 0;

  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum(abs(xv - xx.col(j)) / (xv + xx.col(j)));
    }
  }
  return f;
}

NumericVector wave_hedges_vec(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  NumericVector f(proper_size(nrw, ncl));
  colvec xv(nrw);
  int i, j,k=0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      f[k] = sum(abs(xv - xx.col(j)) / max_elems(xv , xx.col(j)));
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericVector haversine_vec(NumericMatrix x)
{
  const int nrw = x.nrow();
  const int nrw_1 = nrw - 1;
  colvec x0(x.begin(), nrw, false), x1(x.begin() + nrw, nrw, false);
  NumericVector f(proper_size(nrw, nrw));
  colvec ff(f.begin(), f.size(), false);
  colvec ind_col(nrw_1);
  colvec a(nrw_1);
  int i;
  int s = 0, e = 0;
  for (i = 0; i < nrw_1; ++i)
  {
    span ind(i + 1, nrw_1);
    ind_col = x0(ind);
    a = square(sin(0.5 * (x0[i] - ind_col))) + cos(x0[i]) * (cos(ind_col) % square(sin(0.5 * (x1[i] - x1(ind)))));
    a = 2 * asin(sqrt(a));
    e += a.n_elem;
    ff(span(s, e - 1)) = a;
    s += a.n_elem;
  }
  return f;
}

//[[Rcpp::export]]
NumericVector dist_vec(NumericMatrix x, const string method, const bool sqr, const int p)
{
  if (method == "euclidean" || p == 2)
  {
    return euclidean_vec(x, sqr);
  }
  else if (method == "manhattan" || p == 1)
  {
    return manhattan_vec(x);
  }
  else if (method == "maximum")
  {
    return max_vec(x);
  }
  else if (method == "minimum")
  {
    return min_vec(x);
  }
  else if (method == "canberra")
  {
    return canberra_vec(x);
  }
  else if (method == "minkowski")
  {
    return minkowski_vec(x, p);
  }
  else if (method == "bhattacharyya")
  {
    return bhattacharyya_vec(x);
  }
  else if (method == "hellinger")
  {
    return hellinger_vec(x, sqr);
  }
  else if (method == "total_variation")
  {
    return total_variation_vec(x);
  }
  else if (method == "kullback_leibler")
  {
    return kullback_leibler_vec(x);
  }
  else if (method == "jensen_shannon")
  {
    return jensen_shannon_vec(x);
  }
  else if (method == "itakura_saito")
  {
    return itakura_saito_vec(x);
  }
  else if (method == "haversine")
  {
    return haversine_vec(x);
  }
  else if (method == "chi_square")
  {
    return chi_square_vec(x);
  }
  else if (method == "sorensen")
  {
    return sorensen_vec(x);
  }
  else if (method == "soergel")
  {
    return soergel_vec(x);
  }
  else if (method == "cosine")
  {
    return cosine_vec(x);
  }
  else if (method == "wave_hedges")
  {
    return wave_hedges_vec(x);
  }
  stop("Unsupported Method: %s", method);
}

RcppExport SEXP Rfast_dist_vec(SEXP xSEXP, SEXP methodSEXP, SEXP sqrSEXP, SEXP pSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<const string>::type method(methodSEXP);
  traits::input_parameter<const bool>::type sqr(sqrSEXP);
  traits::input_parameter<const int>::type p(pSEXP);
  __result = dist_vec(x, method, sqr, p);
  return __result;
  END_RCPP
}