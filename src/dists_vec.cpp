
// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <string>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;

using std::string;

namespace DistVector
{

  static int proper_size(int nrw, int ncl)
  {
    return ncl * (ncl - 1) * 0.5;
  }

  NumericVector euclidean(NumericMatrix x, const bool sqr)
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

  NumericVector manhattan(NumericMatrix x)
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

  NumericVector hellinger(NumericMatrix x, const bool sqr)
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

  NumericVector max(NumericMatrix x)
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

  NumericVector min(NumericMatrix x)
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
        f[k] = tmp[tmp.index_min()];
      }
    }
    return f;
  }

  NumericVector minkowski(NumericMatrix x, const double p)
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

  NumericVector canberra(NumericMatrix x)
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

  NumericVector total_variation(NumericMatrix x)
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
  NumericVector kullback_leibler(NumericMatrix x)
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
  NumericVector bhattacharyya(NumericMatrix x)
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
        f[k] = -log(sum(sqrt(xv % xx.col(j))));
      }
    }
    return f;
  }

  //[[Rcpp::export]]
  NumericVector jeffries_matusita(NumericMatrix x)
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
        f[k] = sqrt(2.0 - 2.0 * sum(sqrt(xv % xx.col(j))));
      }
    }
    return f;
  }

  //[[Rcpp::export]]
  NumericVector itakura_saito(NumericMatrix x)
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
  NumericVector jensen_shannon(NumericMatrix x)
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

  NumericVector cosine(NumericMatrix x)
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

  NumericVector soergel(NumericMatrix x)
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

  NumericVector chi_square(NumericMatrix x)
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

  NumericVector sorensen(NumericMatrix x)
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

  NumericVector harmonic_mean(NumericMatrix x)
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
        f[k] = 2.0 * dot(xv, xx.col(j)) / sum(xv + xx.col(j));
      }
    }
    return f;
  }

  NumericVector wave_hedges(NumericMatrix x)
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
        f[k] = sum(abs(xv - xx.col(j)) / max_elems(xv, xx.col(j)));
      }
    }
    return f;
  }

  NumericVector motyka(NumericMatrix x)
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
        f[k] = 1.0 - sum_min_elems(xv, xx.col(j)) / sum(xv + xx.col(j));
      }
    }
    return f;
  }

  //[[Rcpp::export]]
  NumericVector haversine(NumericMatrix x)
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
}

//[[Rcpp::export]]
NumericVector dist_vec(NumericMatrix x, const string method, const bool sqr, const int p)
{
  if (method == "euclidean" || p == 2)
  {
    return DistVector::euclidean(x, sqr);
  }
  else if (method == "manhattan" || p == 1)
  {
    return DistVector::manhattan(x);
  }
  else if (method == "maximum")
  {
    return DistVector::max(x);
  }
  else if (method == "minimum")
  {
    return DistVector::min(x);
  }
  else if (method == "canberra")
  {
    return DistVector::canberra(x);
  }
  else if (method == "minkowski")
  {
    return DistVector::minkowski(x, p);
  }
  else if (method == "bhattacharyya")
  {
    return DistVector::bhattacharyya(x);
  }
  else if (method == "hellinger")
  {
    return DistVector::hellinger(x, sqr);
  }
  else if (method == "total_variation")
  {
    return DistVector::total_variation(x);
  }
  else if (method == "kullback_leibler")
  {
    return DistVector::kullback_leibler(x);
  }
  else if (method == "jensen_shannon")
  {
    return DistVector::jensen_shannon(x);
  }
  else if (method == "itakura_saito")
  {
    return DistVector::itakura_saito(x);
  }
  else if (method == "haversine")
  {
    return DistVector::haversine(x);
  }
  else if (method == "chi_square")
  {
    return DistVector::chi_square(x);
  }
  else if (method == "sorensen")
  {
    return DistVector::sorensen(x);
  }
  else if (method == "soergel")
  {
    return DistVector::soergel(x);
  }
  else if (method == "cosine")
  {
    return DistVector::cosine(x);
  }
  else if (method == "wave_hedges")
  {
    return DistVector::wave_hedges(x);
  }
  else if (method == "motyka")
  {
    return DistVector::motyka(x);
  }
  else if (method == "harmonic_mean")
  {
    return DistVector::harmonic_mean(x);
  }
  else if (method == "jeffries_matusita")
  {
    return DistVector::jeffries_matusita(x);
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