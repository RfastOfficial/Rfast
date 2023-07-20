
// Author: Manos Papadakis

//[[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;

double total_euclidean_dist(NumericMatrix x, const bool sqr)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a = 0;
  int i, j;
  if (sqr)
  {
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(square(xx.col(j) - xv));
      }
    }
  }
  else
  {
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += std::sqrt(sum(square(xv - xx.col(j))));
      }
    }
  }
  return a;
}

double total_manhattan_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a = 0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += sum(abs(xv - xx.col(j)));
    }
  }
  return a;
}

double total_hellinger_dist(NumericMatrix x, const bool sqr)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  const double p = 1.0 / std::sqrt(2.0);
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a = 0;
  int i, j;
  if (sqr)
  {
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(square(xv - xx.col(j))) * 0.5;
      }
    }
  }
  else
  {
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += p * std::sqrt(sum(square(xv - xx.col(j))));
      }
    }
  }
  return a;
}

double total_max_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw), tmp(nrw);
  double a = 0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      tmp = abs(xv - xx.col(j));
      a += tmp.at(tmp.index_max());
    }
  }
  return a;
}

double total_min_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a = 0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      xv = abs(xx.col(j) - xv);
      a += xv.at(xv.index_min());
    }
  }
  return a;
}

double total_minkowski_dist(NumericMatrix x, const double p)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  const double p_1 = 1.0 / p;
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a = 0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += pow(sum_with<std::pow, colvec>(abs(xv - xx.col(j)), p), p_1);
    }
  }
  return a;
}

double total_canberra_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw), absx(nrw);
  mat x_abs = abs(x);
  double a = 0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    absx = x_abs.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += sum(abs(xv - xx.col(j)) / (absx + x_abs.col(j)));
    }
  }
  return a;
}

double total_total_variation_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a = 0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += 0.5 * sum(abs(xv - xx.col(j)));
    }
  }
  return a;
}

double total_kullback_leibler_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  NumericMatrix log_x(nrw, ncl);
  mat xx(x.begin(), nrw, ncl, false), log_xx(log_x.begin(), nrw, ncl, false);
  colvec xv(nrw), log_xv(nrw);
  double a = 0;
  int i, j;
  fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    log_xv = log_xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += sum((xv - xx.col(j)) % (log_xv - log_xx.col(j)));
    }
  }
  return a;
}

double total_bhattacharyya_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a = 0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += sum(sqrt(abs(xv % xx.col(j))));
    }
  }
  return a;
}

double total_itakura_saito_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  NumericMatrix log_x(nrw, ncl);
  mat xx(x.begin(), nrw, ncl, false), log_xx(log_x.begin(), nrw, ncl, false);
  colvec xv(nrw), log_xv(nrw);
  double a = 0;
  int i, j;
  fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    log_xv = log_xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += sum((xv - xx.col(j) - log_xv - log_xx.col(j)) - 1);
    }
  }
  return a;
}

//[[Rcpp::export]]
double total_haversine_dist(NumericMatrix x)
{
  const int nrw = x.nrow();
  const int nrw_1 = nrw - 1;
  colvec x0(x.begin(), nrw, false), x1(x.begin() + nrw, nrw, false);
  colvec ind_col(nrw_1);
  double a = 0;
  int i;

  for (i = 0; i < nrw_1; ++i)
  {
    span ind(i + 1, nrw_1);
    ind_col = x0(ind);
    a += accu(2 * asin(sqrt(square(sin(0.5 * (x0[i] - ind_col))) + cos(x0[i]) * (cos(ind_col) % square(sin(0.5 * (x1[i] - x1(ind))))))));
  }
  return a;
}

double total_jensen_shannon_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false), log_xx(log_x.begin(), nrw, ncl, false);
  colvec xv(nrw), log_xv(nrw);
  double a;
  constexpr double log2 = std::log(2);
  int i, j;
  fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());

  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    log_xv = log_xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += sum_with_condition<double, check_if_is_finite, colvec>((xv + xx.col(j)) % (log2 - arma::log(xv + xx.col(j))) + xv % log_xv + xx.col(j) % log_xx.col(j));
    }
  }
  return a;
}
double total_cosine_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw), norm_x = euclidean_norm(xx);
  int i, j, k = 0;

  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    double normx = norm_x[i];
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      a += dot(xv, xx.col(j)) / (normx * norm_x[j]);
    }
  }
  return a;
}

double total_soergel_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      a += sum(abs(xv - xx.col(j))) / sum_max_elems(xv, xx.col(j));
    }
  }
  return a;
}

double total_chi_square_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a;
  int i, j, k = 0;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      a += sum(square(xv - xx.col(j)) / (xv + xx.col(j)));
    }
  }
  return a;
}

double total_sorensen_dist(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a;
  int i, j, k = 0;

  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j, ++k)
    {
      a += sum(abs(xv - xx.col(j)) / (xv + xx.col(j)));
    }
  }
  return a;
}

double total_dists(NumericMatrix x, const string method, const bool sqr, const int p)
{
  if (method == "euclidean" || p == 2)
  {
    return total_euclidean_dist(x, sqr);
  }
  else if (method == "manhattan" || p == 1)
  {
    return total_manhattan_dist(x);
  }
  else if (method == "maximum")
  {
    return total_max_dist(x);
  }
  else if (method == "minimum")
  {
    return total_min_dist(x);
  }
  else if (method == "canberra")
  {
    return total_canberra_dist(x);
  }
  else if (method == "minkowski")
  {
    return total_minkowski_dist(x, p);
  }
  else if (method == "bhattacharyya")
  {
    return total_bhattacharyya_dist(x);
  }
  else if (method == "hellinger")
  {
    return total_hellinger_dist(x, sqr);
  }
  else if (method == "total_variation")
  {
    return total_total_variation_dist(x);
  }
  else if (method == "kullback_leibler")
  {
    return total_kullback_leibler_dist(x);
  }
  else if (method == "jensen_shannon")
  {
    return total_jensen_shannon_dist_vec(x);
  }
  else if (method == "itakura_saito")
  {
    return total_itakura_saito_dist(x);
  }
  else if (method == "haversine")
  {
    return total_haversine_dist(x);
  }
  else if (method == "chi_square")
  {
    return total_chi_square_dist(x);
  }
  else if (method == "sorensen")
  {
    return total_sorensen_dist(x);
  }
  else if (method == "soergel")
  {
    return total_soergel_dist(x);
  }
  else if (method == "cosine")
  {
    return total_cosine_dist(x);
  }
  stop("Unsupported Method: %s", method);
}

RcppExport SEXP Rfast_total_dists(SEXP xSEXP, SEXP methodSEXP, SEXP sqrSEXP, SEXP pSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<const string>::type method(methodSEXP);
  traits::input_parameter<const bool>::type sqr(sqrSEXP);
  traits::input_parameter<const int>::type p(pSEXP);
  __result = total_dists(x, method, sqr, p);
  return __result;
  END_RCPP
}