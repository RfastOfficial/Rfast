
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

double total_manhattan(NumericMatrix x)
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

double total_hellinger(NumericMatrix x, const bool sqr)
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

double total_max(NumericMatrix x)
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

double total_min(NumericMatrix x)
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

double total_minkowski(NumericMatrix x, const double p)
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

double total_canberra(NumericMatrix x)
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

double total_total_variation(NumericMatrix x)
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

double total_kullback_leibler(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
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

double total_bhattacharyya(NumericMatrix x)
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
      a += -log(sum(sqrt(xv % xx.col(j))));
    }
  }
  return a;
}

double total_itakura_saito(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
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
double total_haversine(NumericMatrix x)
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

double total_jensen_shannon(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
  colvec xv(nrw), log_xv(nrw);
  double a=0.0;
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
double total_cosine(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw), norm_x = euclidean_norm(xx);
  int i, j, k = 0;
  double a=0.0;

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

double total_wave_hedges(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a=0.0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += sum(abs(xv - xx.col(j)) / max_elems(xv , xx.col(j)));
    }
  }
  return a;
}

double total_motyka(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a=0.0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += 1.0 - sum_min_elems(xv , xx.col(j)) / sum(xv + xx.col(j));
    }
  }
  return a;
}

double total_harmonic_mean(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a=0.0;
  int i, j;
  for (i = 0; i < ncl - 1; ++i)
  {
    xv = xx.col(i);
    for (j = i + 1; j < ncl; ++j)
    {
      a += 2.0 * dot(xv,xx.col(j)) / sum(xv + xx.col(j));
    }
  }
  return a;
}

double total_soergel(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a=0.0;
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

double total_chi_square(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a=0.0;
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

double total_sorensen(NumericMatrix x)
{
  const int ncl = x.ncol(), nrw = x.nrow();
  mat xx(x.begin(), nrw, ncl, false);
  colvec xv(nrw);
  double a=0.0;
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
    return total_manhattan(x);
  }
  else if (method == "maximum")
  {
    return total_max(x);
  }
  else if (method == "minimum")
  {
    return total_min(x);
  }
  else if (method == "canberra")
  {
    return total_canberra(x);
  }
  else if (method == "minkowski")
  {
    return total_minkowski(x, p);
  }
  else if (method == "bhattacharyya")
  {
    return total_bhattacharyya(x);
  }
  else if (method == "hellinger")
  {
    return total_hellinger(x, sqr);
  }
  else if (method == "total_variation")
  {
    return total_total_variation(x);
  }
  else if (method == "kullback_leibler")
  {
    return total_kullback_leibler(x);
  }
  else if (method == "jensen_shannon")
  {
    return total_jensen_shannon(x);
  }
  else if (method == "itakura_saito")
  {
    return total_itakura_saito(x);
  }
  else if (method == "haversine")
  {
    return total_haversine(x);
  }
  else if (method == "chi_square")
  {
    return total_chi_square(x);
  }
  else if (method == "sorensen")
  {
    return total_sorensen(x);
  }
  else if (method == "soergel")
  {
    return total_soergel(x);
  }
  else if (method == "cosine")
  {
    return total_cosine(x);
  }
  else if (method == "wave_hedges")
  {
    return total_wave_hedges(x);
  }
  else if (method == "motyka")
  {
    return total_motyka(x);
  }
  else if (method == "harmonic_mean")
  {
    return total_harmonic_mean(x);
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