
// Author: Manos Papadakis

//[[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"
#include "Coeff.h"
#include "Dist.h"

using namespace Rcpp;
using namespace arma;

namespace DistTotal
{
	template <class Function>
	double dist_h(NumericMatrix x, Function func)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		colvec xv(nrw);
		double a=0.0;
		size_t i, j;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a += accu(func(xv, xx.col(j)));
			}
		}
		return a;
	}
	template <class Function>
	double dist_h(NumericMatrix x, const double p, Function func)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		NumericMatrix f(ncl, ncl);
		colvec xv(nrw);
		double a=0.0;
		size_t i, j;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a += accu(func(xv, xx.col(j), p));
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
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		colvec xv(nrw);
		double a = 0;
		size_t i, j;
		for (i = 0; i < ncl - 1; ++i)
		{
			xv = sqrt_xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a += -log(Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
			}
		}
		return -a;
	}

	double jeffries_matusita(NumericMatrix x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		colvec xv(nrw);
		double a = 0.0;
		size_t i, j;

		for (i = 0; i < ncl - 1; ++i)
		{
			xv = sqrt_xx.col(i);
			for (j = i + 1; j < ncl; ++j)
			{
				a += sqrt(2.0 - 2.0 * Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
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
				a += sum_with_condition<double, std::isfinite, colvec>(xv / xx.col(j) - (log_xv - log_xx.col(j)) - 1);
			}
		}
		return a;
	}

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

}

double total_dist(NumericMatrix x, const string method, const bool sqr, const int p)
{
	if (method == "canberra")
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
	else if (method == "cosine")
	{
		return DistTotal::cosine(x);
	}
	else if (method == "jeffries_matusita")
	{
		return DistTotal::jeffries_matusita(x);
	}
	else if (method == "manhattan" || p == 2)
	{
		return DistTotal::dist_h(x, Dist::manhattan);
	}
	else if (method == "chi_square")
	{
		return DistTotal::dist_h(x, Dist::chi_square);
	}
	else if (method == "soergel")
	{
		return DistTotal::dist_h(x, Dist::soergel);
	}
	else if (method == "kulczynski")
	{
		return DistTotal::dist_h(x, Dist::kulczynski);
	}
	else if (method == "wave_hedges")
	{
		return DistTotal::dist_h(x, Dist::wave_hedges);
	}
	else if (method == "motyka")
	{
		return DistTotal::dist_h(x, Dist::motyka);
	}
	else if (method == "harmonic_mean")
	{
		return DistTotal::dist_h(x, Dist::harmonic_mean);
	}
	else if (method == "total_variation")
	{
		return DistTotal::dist_h(x, Dist::total_variation);
	}
	else if (method == "sorensen")
	{
		return DistTotal::dist_h(x, Dist::sorensen);
	}
	else if (method == "euclidean" || p == 1)
	{
		return sqr ? DistTotal::dist_h(x, Dist::euclidean<true>) : DistTotal::dist_h(x, Dist::euclidean<false>);
	}
	else if (method == "maximum")
	{
		return DistTotal::dist_h(x, Dist::max);
	}
	else if (method == "minimum")
	{
		return DistTotal::dist_h(x, Dist::min);
	}
	else if (method == "hellinger")
	{
		return sqr ? DistTotal::dist_h(x, 0.5, Dist::hellinger<true>) : DistTotal::dist_h(x, 1.0 / std::sqrt(2.0), Dist::hellinger<false>);
	}
	else if (method == "gower")
	{
		return DistTotal::dist_h(x, 1.0 / x.nrow(), Dist::gower);
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
	__result = total_dist(x, method, sqr, p);
	return __result;
	END_RCPP
}