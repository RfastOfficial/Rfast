
// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"
#include "Coeff.h"
#include "Dist.h"
#include <string>

using namespace Rcpp;
using namespace arma;
using std::string;

namespace Dist
{
	template <class Function>
	NumericMatrix dist_h(NumericMatrix x, Function func)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericMatrix f(ncl, ncl);
		colvec xv(nrw);
		double a;
		size_t i, j;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a = func(xv, xx.col(j));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}
	template <class Function>
	NumericMatrix dist_h(NumericMatrix x, const double p, Function func)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericMatrix f(ncl, ncl);
		colvec xv(nrw);
		double a;
		size_t i, j;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a = func(xv, xx.col(j), p);
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix minkowski(NumericMatrix x, const double p)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericMatrix f(ncl, ncl);
		colvec xv(nrw);
		double a;
		size_t i, j;
		const double p_1 = 1.0 / p;

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a = pow(sum_with<std::pow, colvec>(abs(xv - xx.col(j)), p), p_1);
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix canberra(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericMatrix f(ncl, ncl);
		colvec xv(nrw), absx(nrw);
		mat x_abs = abs(xx);
		double a;
		size_t i, j;

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			absx = x_abs.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a = sum(abs(xv - xx.col(j)) / (absx + x_abs.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix cosine(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericMatrix f(ncl, ncl);
		colvec xv(nrw), norm_x = euclidean_norm(xx).t();
		double a;
		size_t i, j;

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			double normx = norm_x[i];
			for (j = i + 1; j < ncl; ++j)
			{
				a = dot(xv, xx.col(j)) / (normx * norm_x[j]);
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix kullback_leibler(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl), log_x(nrw, ncl);
		mat xx(x.begin(), nrw, ncl, false), log_xx(log_x.begin(), nrw, ncl, false);
		colvec xv(nrw), log_xv(nrw);
		double a;
		size_t i, j;

		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			log_xv = log_xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a = sum_with_condition<double, std::isfinite, colvec>((xv - xx.col(j)) % (log_xv - log_xx.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix jensen_shannon(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl);
		mat xx(x.begin(), nrw, ncl, false);
		colvec xv(nrw), xlogx_xv(nrw);
		mat xlogx = xx % arma::log(xx);
		double a;
		const double log0_5 = std::log(0.5);
		size_t i, j;

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			xlogx_xv = xlogx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{

				a = sum_with_condition<double, check_if_is_finite, colvec>(xlogx_xv + xlogx_xv.col(j) - (arma::log(xv + xx.col(j)) + log0_5) % (xv + xx.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix bhattacharyya(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericMatrix f(ncl, ncl);
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		colvec xv(nrw);
		double a;
		size_t i, j;

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = sqrt_xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a = -log(Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix jeffries_matusita(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericMatrix f(ncl, ncl);
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		colvec xv(nrw);
		double a;
		size_t i, j;

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = sqrt_xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a = sqrt(2.0 - 2.0 * Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix itakura_saito(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl), log_x(nrw, ncl);
		mat xx(x.begin(), nrw, ncl, false), log_xx(log_x.begin(), nrw, ncl, false);
		colvec xv(nrw), log_xv(nrw);
		double a;
		size_t i, j;
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			log_xv = log_xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a = sum_with_condition<double, std::isfinite, colvec>(xv / xx.col(j) - (log_xv - log_xx.col(j)) - 1);
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix haversine(NumericMatrix x)
	{
		const int nrw = x.nrow();
		const int nrw_1 = nrw - 1;
		colvec x0(x.begin(), nrw, false), x1(x.begin() + nrw, nrw, false);
		NumericMatrix f(nrw, nrw);
		mat ff(f.begin(), nrw, nrw, false);
		colvec ind_col(nrw_1);
		colvec a(nrw_1);
		int i;

		for (i = 0; i < nrw_1; ++i)
		{
			span ind(i + 1, nrw_1);
			ind_col = x0(ind);
			a = square(sin(0.5 * (x0[i] - ind_col))) + cos(x0[i]) * (cos(ind_col) % square(sin(0.5 * (x1[i] - x1(ind)))));
			a = 2 * asin(sqrt(a));
			ff(i, ind) = a.t();
			ff(ind, i) = a;
		}
		return f;
	}

}

NumericMatrix dist(NumericMatrix x, const string method, const bool sqr, const int p)
{
	if (method == "canberra")
	{
		return Dist::canberra(x);
	}
	else if (method == "minkowski")
	{
		return Dist::minkowski(x, p);
	}
	else if (method == "bhattacharyya")
	{
		return Dist::bhattacharyya(x);
	}
	else if (method == "kullback_leibler")
	{
		return Dist::kullback_leibler(x);
	}
	else if (method == "jensen_shannon")
	{
		return Dist::jensen_shannon(x);
	}
	else if (method == "itakura_saito")
	{
		return Dist::itakura_saito(x);
	}
	else if (method == "haversine")
	{
		return Dist::haversine(x);
	}
	else if (method == "cosine")
	{
		return Dist::cosine(x);
	}
	else if (method == "jeffries_matusita")
	{
		return Dist::jeffries_matusita(x);
	}
	else if (method == "manhattan" || p == 2)
	{
		return Dist::dist_h(x, Dist::manhattan);
	}
	else if (method == "chi_square")
	{
		return Dist::dist_h(x, Dist::chi_square);
	}
	else if (method == "soergel")
	{
		return Dist::dist_h(x, Dist::soergel);
	}
	else if (method == "kulczynski")
	{
		return Dist::dist_h(x, Dist::kulczynski);
	}
	else if (method == "wave_hedges")
	{
		return Dist::dist_h(x, Dist::wave_hedges);
	}
	else if (method == "motyka")
	{
		return Dist::dist_h(x, Dist::motyka);
	}
	else if (method == "harmonic_mean")
	{
		return Dist::dist_h(x, Dist::harmonic_mean);
	}
	else if (method == "total_variation")
	{
		return Dist::dist_h(x, Dist::total_variation);
	}
	else if (method == "sorensen")
	{
		return Dist::dist_h(x, Dist::sorensen);
	}
	else if (method == "euclidean" || p == 1)
	{
		return sqr ? Dist::dist_h(x, Dist::euclidean<true>) : Dist::dist_h(x, Dist::euclidean<false>);
	}
	else if (method == "maximum")
	{
		return Dist::dist_h(x, Dist::max);
	}
	else if (method == "minimum")
	{
		return Dist::dist_h(x, Dist::min);
	}
	else if (method == "hellinger")
	{
		return sqr ? Dist::dist_h(x, 0.5, Dist::hellinger<true>) : Dist::dist_h(x, 1.0 / std::sqrt(2.0), Dist::hellinger<false>);
	}
	else if (method == "gower")
	{
		return Dist::dist_h(x, 1.0 / x.nrow(), Dist::gower);
	}
	stop("Unsupported Method: %s", method);
}

RcppExport SEXP Rfast_dist(SEXP xSEXP, SEXP methodSEXP, SEXP sqrSEXP, SEXP pSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type sqr(sqrSEXP);
	traits::input_parameter<const int>::type p(pSEXP);
	__result = dist(x, method, sqr, p);
	return __result;
	END_RCPP
}

namespace DistVector
{

	template <class Function>
	NumericVector dist_h(NumericMatrix x, Function func)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericVector f(proper_size(nrw, ncl));
		colvec xv(nrw);
		size_t i, j, k = 0;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			for (j = i + 1; j < ncl; ++j, ++k)
			{
				f[k] = func(xv, xx.col(j));
			}
		}
		return f;
	}

	template <class Function>
	NumericVector dist_h(NumericMatrix x, const double p, Function func)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericVector f(proper_size(nrw, ncl));
		colvec xv(nrw);
		size_t i, j, k = 0;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			for (j = i + 1; j < ncl; ++j, ++k)
			{
				f[k] = func(xv, xx.col(j), p);
			}
		}
		return f;
	}

	NumericVector minkowski(NumericMatrix x, const double p)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		const double p_1 = 1.0 / p;
		mat xx(x.begin(), nrw, ncl, false);
		NumericVector f(proper_size(nrw, ncl));
		colvec xv(nrw);
		size_t i, j, k = 0;
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
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericVector f(proper_size(nrw, ncl));
		colvec xv(nrw), absx(nrw);
		mat x_abs = abs(x);
		size_t i, j, k = 0;
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

	NumericVector kullback_leibler(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
		colvec xv(nrw), log_xv(nrw);
		size_t i, j, k = 0;
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

	NumericVector bhattacharyya(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericVector f(proper_size(nrw, ncl));
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		colvec xv(nrw);
		size_t i, j, k = 0;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = sqrt_xx.col(i);
			for (j = i + 1; j < ncl; ++j, ++k)
			{
				f[k] = -log(Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
			}
		}
		return f;
	}

	NumericVector jeffries_matusita(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericVector f(proper_size(nrw, ncl));
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		colvec xv(nrw);
		size_t i, j, k = 0;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = sqrt_xx.col(i);
			for (j = i + 1; j < ncl; ++j, ++k)
			{
				f[k] = sqrt(2.0 - 2.0 * Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
			}
		}
		return f;
	}

	NumericVector itakura_saito(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
		colvec xv(nrw), log_xv(nrw);
		size_t i, j, k = 0;
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

	NumericVector jensen_shannon(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
		mat xlogx = xx % arma::log(xx);
		colvec xv(nrw), xlogx_xv(nrw);
		const double log0_5 = std::log(0.5);
		size_t i, j, k = 0;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			xlogx_xv = xlogx.col(i);
			for (j = i + 1; j < ncl; ++j, ++k)
			{
				f[k] = sum_with_condition<double, check_if_is_finite, colvec>(xlogx_xv + xlogx_xv.col(j) - (arma::log(xv + xx.col(j)) + log0_5) % (xv + xx.col(j)));
			}
		}
		return f;
	}

	NumericVector cosine(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false);
		colvec xv(nrw), norm_x = euclidean_norm(xx).t();
		size_t i, j, k = 0;

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
	if (method == "canberra")
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
	else if (method == "cosine")
	{
		return DistVector::cosine(x);
	}
	else if (method == "jeffries_matusita")
	{
		return DistVector::jeffries_matusita(x);
	}
	else if (method == "manhattan" || p == 2)
	{
		return DistVector::dist_h(x, Dist::manhattan);
	}
	else if (method == "chi_square")
	{
		return DistVector::dist_h(x, Dist::chi_square);
	}
	else if (method == "soergel")
	{
		return DistVector::dist_h(x, Dist::soergel);
	}
	else if (method == "kulczynski")
	{
		return DistVector::dist_h(x, Dist::kulczynski);
	}
	else if (method == "wave_hedges")
	{
		return DistVector::dist_h(x, Dist::wave_hedges);
	}
	else if (method == "motyka")
	{
		return DistVector::dist_h(x, Dist::motyka);
	}
	else if (method == "harmonic_mean")
	{
		return DistVector::dist_h(x, Dist::harmonic_mean);
	}
	else if (method == "total_variation")
	{
		return DistVector::dist_h(x, Dist::total_variation);
	}
	else if (method == "sorensen")
	{
		return DistVector::dist_h(x, Dist::sorensen);
	}
	else if (method == "euclidean" || p == 1)
	{
		return sqr ? DistVector::dist_h(x, Dist::euclidean<true>) : DistVector::dist_h(x, Dist::euclidean<false>);
	}
	else if (method == "maximum")
	{
		return DistVector::dist_h(x, Dist::max);
	}
	else if (method == "minimum")
	{
		return DistVector::dist_h(x, Dist::min);
	}
	else if (method == "hellinger")
	{
		return sqr ? DistVector::dist_h(x, 0.5, Dist::hellinger<true>) : DistVector::dist_h(x, 1.0 / std::sqrt(2.0), Dist::hellinger<false>);
	}
	else if (method == "gower")
	{
		return DistVector::dist_h(x, 1.0 / x.nrow(), Dist::gower);
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