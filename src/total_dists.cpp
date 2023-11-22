
// Author: Manos Papadakis

//[[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"
#include "Coeff.h"

using namespace Rcpp;
using namespace arma;

namespace DistTotal
{

  double euclidean(NumericMatrix x, const bool sqr)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0;
    size_t i, j;
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

  double manhattan(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0;
    size_t i, j;
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

  double hellinger(NumericMatrix x, const bool sqr)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    const double p = 1.0 / std::sqrt(2.0);
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0;
    size_t i, j;
    if (sqr)
    {
      for (i = 0; i < ncl - 1; ++i)
      {
        xv = xx.col(i);
        for (j = i + 1; j < ncl; ++j)
        {
          a += sum(square(xv - xx.col(j)));
        }
      }
      a *= 0.5;
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
      a *= p;
    }
    return a;
  }

  double max(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw), tmp(nrw);
    double a = 0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        tmp = abs(xv - xx.col(j));
        a += tmp[tmp.index_max()];
      }
    }
    return a;
  }

  double min(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw), tmp(nrw);
    double a = 0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        tmp = abs(xx.col(j) - xv);
        a += tmp[tmp.index_min()];
      }
    }
    return a;
  }

  double minkowski(NumericMatrix x, const double p)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    const double p_1 = 1.0 / p;
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0;
    size_t i, j;
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

  double canberra(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw), absx(nrw);
    mat x_abs = abs(x);
    double a = 0;
    size_t i, j;
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

  //[[Rcpp::export]]
  double gower(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    const double p = 1.0 / nrw;
    NumericMatrix f(ncl, ncl);
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw), log_xv(nrw);
    double a = 0.0;
    size_t i, j;

    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(abs(xv - xx.col(j)));
      }
    }
    return a * p;
  }

  double kulczynski(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0.0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(abs(xv - xx.col(j))) / sum_min_elems(xv, xx.col(j));
      }
    }
    return a;
  }

  double total_variation(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(abs(xv - xx.col(j)));
      }
    }
    return a * 0.5;
  }

  double kullback_leibler(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
    colvec xv(nrw), log_xv(nrw);
    double a = 0;
    size_t i, j;
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

  double bhattacharyya(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += -log(Coeff::bhattacharyya(xv, xx.col(j)));
      }
    }
    return -a;
  }

  //[[Rcpp::export]]
  double jeffries_matusita(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0.0;
    size_t i, j;

    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sqrt(2.0 - 2.0 * Coeff::bhattacharyya(xv, xx.col(j)));
      }
    }
    return a;
  }

  double itakura_saito(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
    colvec xv(nrw), log_xv(nrw);
    double a = 0;
    size_t i, j;
    fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      log_xv = log_xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum_with_condition<double,std::isfinite, colvec>(xv / xx.col(j) - (log_xv - log_xx.col(j)) - 1);
      }
    }
    return a;
  }

  //[[Rcpp::export]]
  double haversine(NumericMatrix x)
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

  double jensen_shannon(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    NumericMatrix f(ncl, ncl);
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw), xlogx_xv(nrw);
    mat xlogx = xx % arma::log(xx);
    double a = 0.0;
    const double log0_5 = std::log(0.5);
    size_t i, j;

    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      xlogx_xv = xlogx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {

        a += sum_with_condition<double, check_if_is_finite, colvec>(xlogx_xv + xlogx_xv.col(j) - (arma::log(xv + xx.col(j)) + log0_5) % (xv + xx.col(j)));
      }
    }
    return a;
  }

  double cosine(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw), norm_x = euclidean_norm(xx);
    size_t i, j;
    double a = 0.0;

    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      double normx = norm_x[i];
      for (j = i + 1; j < ncl; ++j)
      {
        a += dot(xv, xx.col(j)) / (normx * norm_x[j]);
      }
    }
    return a;
  }

  double wave_hedges(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0.0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(abs(xv - xx.col(j)) / max_elems(xv, xx.col(j)));
      }
    }
    return a;
  }

  double motyka(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0.0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += 1.0 - sum_min_elems(xv, xx.col(j)) / sum(xv + xx.col(j));
      }
    }
    return a;
  }

  double harmonic_mean(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0.0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += dot(xv, xx.col(j)) / sum(xv + xx.col(j));
      }
    }
    return a * 2.0;
  }

  double soergel(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0.0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(abs(xv - xx.col(j))) / sum_max_elems(xv, xx.col(j));
      }
    }
    return a;
  }

  double chi_square(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0.0;
    size_t i, j;
    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(square(xv - xx.col(j)) / (xv + xx.col(j)));
      }
    }
    return a;
  }

  double sorensen(NumericMatrix x)
  {
    const size_t ncl = x.ncol(), nrw = x.nrow();
    mat xx(x.begin(), nrw, ncl, false);
    colvec xv(nrw);
    double a = 0.0;
    size_t i, j;

    for (i = 0; i < ncl - 1; ++i)
    {
      xv = xx.col(i);
      for (j = i + 1; j < ncl; ++j)
      {
        a += sum(abs(xv - xx.col(j)) / (xv + xx.col(j)));
      }
    }
    return a;
  }

}

double total_dists(NumericMatrix x, const string method, const bool sqr, const int p)
{
  if (method == "euclidean" || p == 2)
  {
    return DistTotal::euclidean(x, sqr);
  }
  else if (method == "manhattan" || p == 1)
  {
    return DistTotal::manhattan(x);
  }
  else if (method == "maximum")
  {
    return DistTotal::max(x);
  }
  else if (method == "minimum")
  {
    return DistTotal::min(x);
  }
  else if (method == "canberra")
  {
    return DistTotal::canberra(x);
  }
  else if (method == "minkowski")
  {
    return DistTotal::minkowski(x, p);
  }
  else if (method == "bhattacharyya")
  {
    return DistTotal::bhattacharyya(x);
  }
  else if (method == "hellinger")
  {
    return DistTotal::hellinger(x, sqr);
  }
  else if (method == "total_variation")
  {
    return DistTotal::total_variation(x);
  }
  else if (method == "kullback_leibler")
  {
    return DistTotal::kullback_leibler(x);
  }
  else if (method == "jensen_shannon")
  {
    return DistTotal::jensen_shannon(x);
  }
  else if (method == "itakura_saito")
  {
    return DistTotal::itakura_saito(x);
  }
  else if (method == "haversine")
  {
    return DistTotal::haversine(x);
  }
  else if (method == "chi_square")
  {
    return DistTotal::chi_square(x);
  }
  else if (method == "sorensen")
  {
    return DistTotal::sorensen(x);
  }
  else if (method == "soergel")
  {
    return DistTotal::soergel(x);
  }
  else if (method == "cosine")
  {
    return DistTotal::cosine(x);
  }
  else if (method == "wave_hedges")
  {
    return DistTotal::wave_hedges(x);
  }
  else if (method == "motyka")
  {
    return DistTotal::motyka(x);
  }
  else if (method == "harmonic_mean")
  {
    return DistTotal::harmonic_mean(x);
  }
  else if (method == "jeffries_matusita")
  {
    return DistTotal::jeffries_matusita(x);
  }
  else if (method == "gower")
  {
    return DistTotal::gower(x);
  }
  else if (method == "kulczynski")
  {
    return DistTotal::kulczynski(x);
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