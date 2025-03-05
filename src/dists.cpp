
// Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast/Dist.h"
#include "Rfast.h"
#include "Coeff.h"
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;
using std::string;

namespace Dist
{
	template <class Function>
	void dist_inner(mat &xx, colvec &xv, size_t i, size_t ncl, size_t nrw, mat &ff, Function func)
	{
		for (size_t j = i + 1; j < ncl; ++j)
		{
			colvec y(xx.begin_col(j), nrw, false);
			double a = func(xv, y);
			ff(i, j) = a;
			ff(j, i) = a;
		}
	}

	template <class Function>
	NumericMatrix dist_h(NumericMatrix &x, Function func, const bool parallel = false)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl);
		mat xx(x.begin(), nrw, ncl, false), ff(f.begin(), ncl, ncl, false);
		if (parallel)
		{
#ifdef _OPENMP
	#pragma omp parallel for
#endif
			for (size_t i = 0; i < ncl - 1; ++i)
			{
				colvec xv(xx.begin_col(i), nrw, false);
				dist_inner<Function>(xx, xv, i, ncl, nrw, ff, func);
			}
		}
		else
		{
			for (size_t i = 0; i < ncl - 1; ++i)
			{
				colvec xv(xx.begin_col(i), nrw, false);
				dist_inner<Function>(xx, xv, i, ncl, nrw, ff, func);
			}
		}
		return f;
	}

	template <class Function>
	NumericMatrix dist_h(NumericMatrix &x, const double p, Function func, const bool parallel = false)
	{
		auto func2 = [&](colvec &x, colvec &y)
		{ return func(x, y, p); };
		return dist_h(x, func2, parallel);
	}

	NumericMatrix canberra(NumericMatrix &x, const bool parallel = false)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl);
		mat xx(x.begin(), nrw, ncl, false), ff(f.begin(), ncl, ncl, false);
		mat x_abs = abs(xx);

		if (parallel)
		{
#ifdef _OPENMP
	#pragma omp parallel for
#endif
			for (size_t i = 0; i < ncl - 1; ++i)
			{
				colvec xv(xx.begin_col(i), nrw, false);
				colvec absx(x_abs.begin_col(i), nrw, false);
				for (size_t j = i + 1; j < ncl; ++j)
				{
					double a = sum(abs(xv - xx.col(j)) / (absx + x_abs.col(j)));
					ff(i, j) = a;
					ff(j, i) = a;
				}
			}
		}
		else
		{
			for (size_t i = 0; i < ncl - 1; ++i)
			{
				colvec xv(xx.begin_col(i), nrw, false);
				colvec absx(x_abs.begin_col(i), nrw, false);
				for (size_t j = i + 1; j < ncl; ++j)
				{
					double a = sum(abs(xv - xx.col(j)) / (absx + x_abs.col(j)));
					ff(i, j) = a;
					ff(j, i) = a;
				}
			}
		}
		return f;
	}

	NumericMatrix cosine(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl);
		mat xx(x.begin(), nrw, ncl, false), ff(f.begin(), ncl, ncl, false);
		colvec norm_x = euclidean_norm(xx).t();
		double a = 0;

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			double normx = norm_x[i];
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a = 2.0 * (1 - dot(xv, xx.col(j)) / (normx * norm_x[j]));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix kullback_leibler(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl), log_x(nrw, ncl);
		mat xx(x.begin(), nrw, ncl, false), log_xx(log_x.begin(), nrw, ncl, false);
		colvec log_xv(nrw);
		double a = 0;

		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec log_xv(log_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a = sum_with_condition<double, std::isfinite, colvec>((xv - xx.col(j)) % (log_xv - log_xx.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix itakura_saito(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl), log_x(nrw, ncl);
		mat xx(x.begin(), nrw, ncl, false), log_xx(log_x.begin(), nrw, ncl, false);
		colvec log_xv(nrw);
		double a = 0;
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec log_xv(log_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a = sum_with_condition<double, std::isfinite, colvec>(xv / xx.col(j) - (log_xv - log_xx.col(j)) - 1);
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix jensen_shannon(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl);
		mat xx(x.begin(), nrw, ncl, false), ff(f.begin(), ncl, ncl, false);
		mat xlogx = xx % arma::log(xx);
		double a = 0;
		const double log0_5 = std::log(0.5);

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec xlogx_xv(xlogx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a = sum_with_condition<double, check_if_is_finite, colvec>(xlogx_xv + xlogx.col(j) - (arma::log(xv + xx.col(j)) + log0_5) % (xv + xx.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix bhattacharyya(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl);
		mat xx(x.begin(), nrw, ncl, false), ff(f.begin(), ncl, ncl, false);
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		double a = 0;

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(sqrt_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a = -log(Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix jeffries_matusita(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl);
		mat xx(x.begin(), nrw, ncl, false), ff(f.begin(), ncl, ncl, false);
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		double a = 0;

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(sqrt_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a = sqrt(2.0 - 2.0 * Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
				f(i, j) = a;
				f(j, i) = a;
			}
		}
		return f;
	}

	NumericMatrix haversine(NumericMatrix &x)
	{
		const size_t nrw = x.nrow(), nrw_1 = nrw - 1;
		NumericMatrix f(nrw, nrw);
		mat ff(f.begin(), nrw, nrw, false);
		colvec x0(x.begin(), nrw, false), x1(x.begin() + nrw, nrw, false), ind_col(nrw_1), a(nrw_1);

		for (size_t i = 0; i < nrw_1; ++i)
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

NumericMatrix dist(NumericMatrix x, const string method, const bool sqr, const int p, const bool parallel)
{
	if (method == "euclidean" || p == 1)
	{
		return sqr ? Dist::dist_h(x, Rfast::Dist::euclidean<false, colvec>, parallel) : Dist::dist_h(x, Rfast::Dist::euclidean<true, colvec>, parallel);
	}
	else if (method == "manhattan" || p == 2)
	{
		return Dist::dist_h(x, Rfast::Dist::manhattan);
	}
	else if (method == "canberra")
	{
		return Dist::canberra(x, parallel);
	}
	else if (method == "minkowski")
	{
		return Dist::dist_h(x, p, Rfast::Dist::minkowski, parallel);
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
	else if (method == "chi_square")
	{
		return Dist::dist_h(x, Rfast::Dist::chi_square);
	}
	else if (method == "soergel")
	{
		return Dist::dist_h(x, Rfast::Dist::soergel);
	}
	else if (method == "kulczynski")
	{
		return Dist::dist_h(x, Rfast::Dist::kulczynski);
	}
	else if (method == "wave_hedges")
	{
		return Dist::dist_h(x, Rfast::Dist::wave_hedges);
	}
	else if (method == "motyka")
	{
		return Dist::dist_h(x, Rfast::Dist::motyka);
	}
	else if (method == "harmonic_mean")
	{
		return Dist::dist_h(x, Rfast::Dist::harmonic_mean);
	}
	else if (method == "total_variation")
	{
		return Dist::dist_h(x, Rfast::Dist::total_variation);
	}
	else if (method == "sorensen")
	{
		return Dist::dist_h(x, Rfast::Dist::sorensen);
	}
	else if (method == "maximum")
	{
		return Dist::dist_h(x, Rfast::Dist::max);
	}
	else if (method == "minimum")
	{
		return Dist::dist_h(x, Rfast::Dist::min);
	}
	else if (method == "hellinger")
	{
		return sqr ? Dist::dist_h(x, 0.5, Rfast::Dist::hellinger<true>) : Dist::dist_h(x, 1.0 / std::sqrt(2.0), Rfast::Dist::hellinger<false>);
	}
	else if (method == "gower")
	{
		return Dist::dist_h(x, 1.0 / x.nrow(), Rfast::Dist::gower);
	}
	stop("Unsupported Method: %s", method);
}

RcppExport SEXP Rfast_dist(SEXP xSEXP, SEXP methodSEXP, SEXP sqrSEXP, SEXP pSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type sqr(sqrSEXP);
	traits::input_parameter<const int>::type p(pSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = dist(x, method, sqr, p, parallel);
	return __result;
	END_RCPP
}

namespace DistVector
{
	template <class Function, bool parallel>
	void dist_inner(mat &xx, colvec &xv, size_t i, size_t ncl, size_t nrw, colvec &ff, size_t &k, Function func)
	{
		if constexpr (parallel)
		{
#pragma omp parallel
			{
				for (size_t j = i + 1; j < ncl; ++j)
				{
					colvec y(xx.begin_col(j), nrw, false);
					ff[k] = func(xv, y);	
#ifdef _OPENMP
	#pragma omp atomic
#endif
					++k;
				}
			}
		}
		else
		{
			for (size_t j = i + 1; j < ncl; ++j)
			{
				colvec y(xx.begin_col(j), nrw, false);
				ff[k++] = func(xv, y);
			}
		}
	}

	template <class Function>
	NumericVector dist_h(NumericMatrix &x, Function func, const bool parallel = false)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		colvec ff(f.begin(), f.size(), false);
		mat xx(x.begin(), nrw, ncl, false);
		size_t k = 0;
		if (parallel)
		{
#ifdef _OPENMP
	#pragma omp parallel for
#endif
			for (size_t i = 0; i < ncl - 1; ++i)
			{
				colvec xv(xx.begin_col(i), nrw, false);
				dist_inner<Function, true>(xx, xv, i, ncl, nrw, ff, k, func);
			}
		}
		else
		{
			for (size_t i = 0; i < ncl - 1; ++i)
			{
				colvec xv(xx.begin_col(i), nrw, false);
				dist_inner<Function, false>(xx, xv, i, ncl, nrw, ff, k, func);
			}
		}
		return f;
	}

	template <class Function>
	NumericVector dist_h(NumericMatrix &x, const double p, Function func, const bool parallel = false)
	{
		auto func2 = [&](colvec &x, colvec &y)
		{ return func(x, y, p); };
		return dist_h(x, func2, parallel);
	}

	NumericVector canberra(NumericMatrix &x, const bool parallel = false)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false);
		colvec ff(f.begin(), f.size(), false);
		mat x_abs = abs(xx);
		size_t k = 0;
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec absx(x_abs.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j, ++k)
			{
				ff[k] = sum(abs(xv - xx.col(j)) / (absx + x_abs.col(j)));
			}
		}
		return f;
	}

	NumericVector kullback_leibler(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
		colvec log_ff(f.begin(), f.size(), false), ff(f.begin(), f.size(), false);
		size_t k = 0;
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec log_xv(log_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j, ++k)
			{
				ff[k] = sum((xv - xx.col(j)) % (log_xv - log_xx.col(j)));
			}
		}
		return f;
	}

	NumericVector bhattacharyya(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false), sqrt_xx(nrw, ncl, fill::none);
		colvec ff(f.begin(), f.size(), false);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		size_t k = 0;
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(sqrt_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j, ++k)
			{
				ff[k] = -log(Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
			}
		}
		return f;
	}

	NumericVector jeffries_matusita(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false), sqrt_xx(nrw, ncl, fill::none);
		colvec ff(f.begin(), f.size(), false);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		size_t k = 0;
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(sqrt_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j, ++k)
			{
				ff[k] = sqrt(2.0 - 2.0 * Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
			}
		}
		return f;
	}

	NumericVector itakura_saito(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
		colvec log_ff(f.begin(), f.size(), false), ff(f.begin(), f.size(), false);
		size_t k = 0;
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec log_xv(log_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j, ++k)
			{
				ff[k] = sum(xv / xx.col(j) - (log_xv - log_xx.col(j)) - 1);
			}
		}
		return f;
	}

	NumericVector jensen_shannon(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false), xlogx = xx % arma::log(xx);
		colvec ff(f.begin(), f.size(), false);
		const double log0_5 = std::log(0.5);
		size_t k = 0;
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec xlogx_xv(xlogx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j, ++k)
			{
				ff[k] = sum_with_condition<double, check_if_is_finite, colvec>(xlogx_xv + xlogx.col(j) - (arma::log(xv + xx.col(j)) + log0_5) % (xv + xx.col(j)));
			}
		}
		return f;
	}

	NumericVector cosine(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericVector f(proper_size(nrw, ncl));
		mat xx(x.begin(), nrw, ncl, false);
		colvec norm_x = euclidean_norm(xx).t(), ff(f.begin(), f.size(), false);
		size_t k = 0;
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			double normx = norm_x[i];
			for (size_t j = i + 1; j < ncl; ++j, ++k)
			{
				ff[k] = dot(xv, xx.col(j)) / (normx * norm_x[j]);
			}
		}
		return f;
	}

	NumericVector haversine(NumericMatrix &x)
	{
		const size_t nrw = x.nrow(), nrw_1 = nrw - 1;
		NumericVector f(proper_size(nrw, nrw));
		colvec x0(x.begin(), nrw, false), x1(x.begin() + nrw, nrw, false), ff(f.begin(), f.size(), false), ind_col(nrw_1), a(nrw_1);
		int s = 0, e = 0;
		for (size_t i = 0; i < nrw_1; ++i)
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
NumericVector dist_vec(NumericMatrix x, const string method, const bool sqr, const int p, const bool parallel)
{
	if (method == "euclidean" || p == 1)
	{
		return sqr ? DistVector::dist_h(x, Rfast::Dist::euclidean<false, colvec>, parallel) : DistVector::dist_h(x, Rfast::Dist::euclidean<true, colvec>, parallel);
	}
	else if (method == "manhattan" || p == 2)
	{
		return DistVector::dist_h(x, Rfast::Dist::manhattan, parallel);
	}
	else if (method == "canberra")
	{
		return DistVector::canberra(x, parallel);
	}
	else if (method == "minkowski")
	{
		return DistVector::dist_h(x, p, Rfast::Dist::minkowski, parallel);
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
	else if (method == "chi_square")
	{
		return DistVector::dist_h(x, Rfast::Dist::chi_square, parallel);
	}
	else if (method == "soergel")
	{
		return DistVector::dist_h(x, Rfast::Dist::soergel, parallel);
	}
	else if (method == "kulczynski")
	{
		return DistVector::dist_h(x, Rfast::Dist::kulczynski, parallel);
	}
	else if (method == "wave_hedges")
	{
		return DistVector::dist_h(x, Rfast::Dist::wave_hedges, parallel);
	}
	else if (method == "motyka")
	{
		return DistVector::dist_h(x, Rfast::Dist::motyka, parallel);
	}
	else if (method == "harmonic_mean")
	{
		return DistVector::dist_h(x, Rfast::Dist::harmonic_mean, parallel);
	}
	else if (method == "total_variation")
	{
		return DistVector::dist_h(x, Rfast::Dist::total_variation, parallel);
	}
	else if (method == "sorensen")
	{
		return DistVector::dist_h(x, Rfast::Dist::sorensen, parallel);
	}
	else if (method == "maximum")
	{
		return DistVector::dist_h(x, Rfast::Dist::max);
	}
	else if (method == "minimum")
	{
		return DistVector::dist_h(x, Rfast::Dist::min);
	}
	else if (method == "hellinger")
	{
		return sqr ? DistVector::dist_h(x, 0.5, Rfast::Dist::hellinger<true>, parallel) : DistVector::dist_h(x, 1.0 / std::sqrt(2.0), Rfast::Dist::hellinger<false>, parallel);
	}
	else if (method == "gower")
	{
		return DistVector::dist_h(x, 1.0 / x.nrow(), Rfast::Dist::gower, parallel);
	}
	stop("Unsupported Method: %s", method);
}

RcppExport SEXP Rfast_dist_vec(SEXP xSEXP, SEXP methodSEXP, SEXP sqrSEXP, SEXP pSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type sqr(sqrSEXP);
	traits::input_parameter<const int>::type p(pSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = dist_vec(x, method, sqr, p, parallel);
	return __result;
	END_RCPP
}

namespace DistTotal
{
	template <class Function, bool parallel>
	double dist_inner(mat &xx, colvec &xv, size_t i, size_t ncl, size_t nrw, Function func)
	{
		double a = 0.0;
		if constexpr (parallel)
		{
#ifdef _OPENMP
	#pragma omp parallel
	{
#endif
				for (size_t j = i + 1; j < ncl; ++j)
				{
					colvec y(xx.begin_col(j), nrw, false);
					double tmp = func(xv, y);	
#ifdef _OPENMP
	#pragma omp atomic
#endif
					a += tmp;
				}
#ifdef _OPENMP
}
#endif
		}
		else
		{
			for (size_t j = i + 1; j < ncl; ++j)
			{
				colvec y(xx.begin_col(j), nrw, false);
				a += func(xv, y);
			}
		}

		return a;
	}

	template <class Function>
	double dist_h(NumericMatrix &x, Function func, const bool parallel = false)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		double a = 0;
		if (parallel)
		{
#ifdef _OPENMP
	#pragma omp parallel for
#endif
			for (size_t i = 0; i < ncl - 1; ++i)
			{
				colvec xv(xx.begin_col(i), nrw, false);
				double tmp = dist_inner<Function,true>(xx, xv, i, ncl, nrw, func);
				
#ifdef _OPENMP
	#pragma omp atomic
#endif
				a += tmp;
			}
		}
		else
		{
			for (size_t i = 0; i < ncl - 1; ++i)
			{
				colvec xv(xx.begin_col(i), nrw, false);
				a += dist_inner<Function,false>(xx, xv, i, ncl, nrw, func);
			}
		}
		return a;
	}

	template <class Function>
	double dist_h(NumericMatrix &x, const double p, Function func, const bool parallel = false)
	{
		auto func2 = [&](colvec &x, colvec &y)
		{ return func(x, y, p); };
		return dist_h(x, func2, parallel);
	}

	double canberra(NumericMatrix &x, const bool parallel = false)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false), x_abs = abs(xx);
		double a = 0;
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec absx(x_abs.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a += sum(abs(xv - xx.col(j)) / (absx + x_abs.col(j)));
			}
		}
		return a;
	}

	double kullback_leibler(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
		colvec log_xv(nrw);
		double a = 0;
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec log_xv(log_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a += sum((xv - xx.col(j)) % (log_xv - log_xx.col(j)));
			}
		}
		return a;
	}

	double bhattacharyya(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false), sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		double a = 0;
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(sqrt_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a += -log(Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
			}
		}
		return a;
	}

	double jeffries_matusita(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
		colvec xv(nrw);
		double a = 0.0;

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(sqrt_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a += sqrt(2.0 - 2.0 * Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j)));
			}
		}
		return a;
	}

	double itakura_saito(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false), log_xx(nrw, ncl, fill::none);
		colvec log_xv(nrw);
		double a = 0;
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec log_xv(log_xx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a += sum_with_condition<double, std::isfinite, colvec>(xv / xx.col(j) - (log_xv - log_xx.col(j)) - 1);
			}
		}
		return a;
	}

	double haversine(NumericMatrix &x)
	{
		const size_t nrw = x.nrow(), nrw_1 = nrw - 1;
		colvec x0(x.begin(), nrw, false), x1(x.begin() + nrw, nrw, false), ind_col(nrw_1);
		double a = 0;

		for (size_t i = 0; i < nrw_1; ++i)
		{
			span ind(i + 1, nrw_1);
			ind_col = x0(ind);
			a += accu(2 * asin(sqrt(square(sin(0.5 * (x0[i] - ind_col))) + cos(x0[i]) * (cos(ind_col) % square(sin(0.5 * (x1[i] - x1(ind))))))));
		}
		return a;
	}

	double jensen_shannon(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		NumericMatrix f(ncl, ncl);
		mat xx(x.begin(), nrw, ncl, false);
		colvec xlogx_xv(nrw);
		mat xlogx = xx % arma::log(xx);
		double a = 0.0;
		const double log0_5 = std::log(0.5);

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			colvec xlogx_xv(xlogx.begin_col(i), nrw, false);
			for (size_t j = i + 1; j < ncl; ++j)
			{

				a += sum_with_condition<double, check_if_is_finite, colvec>(xlogx_xv + xlogx.col(j) - (arma::log(xv + xx.col(j)) + log0_5) % (xv + xx.col(j)));
			}
		}
		return a;
	}

	double cosine(NumericMatrix &x)
	{
		const size_t ncl = x.ncol(), nrw = x.nrow();
		mat xx(x.begin(), nrw, ncl, false);
		colvec norm_x = euclidean_norm(xx).t();
		double a = 0.0;

		for (size_t i = 0; i < ncl - 1; ++i)
		{
			colvec xv(xx.begin_col(i), nrw, false);
			double normx = norm_x[i];
			for (size_t j = i + 1; j < ncl; ++j)
			{
				a += 2.0 * (1 - dot(xv, xx.col(j)) / (normx * norm_x[j]));
			}
		}
		return a;
	}

}

double total_dist(NumericMatrix x, const string method, const bool sqr, const int p, const bool parallel)
{
	if (method == "euclidean" || p == 1)
	{
		return sqr ? DistTotal::dist_h(x, Rfast::Dist::euclidean<false, colvec>, parallel) : DistTotal::dist_h(x, Rfast::Dist::euclidean<true, colvec>, parallel);
	}
	else if (method == "manhattan" || p == 2)
	{
		return DistTotal::dist_h(x, Rfast::Dist::manhattan, parallel);
	}
	else if (method == "canberra")
	{
		return DistTotal::canberra(x, parallel);
	}
	else if (method == "minkowski")
	{
		return DistTotal::dist_h(x, p, Rfast::Dist::minkowski, parallel);
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
	else if (method == "chi_square")
	{
		return DistTotal::dist_h(x, Rfast::Dist::chi_square, parallel);
	}
	else if (method == "soergel")
	{
		return DistTotal::dist_h(x, Rfast::Dist::soergel, parallel);
	}
	else if (method == "kulczynski")
	{
		return DistTotal::dist_h(x, Rfast::Dist::kulczynski, parallel);
	}
	else if (method == "wave_hedges")
	{
		return DistTotal::dist_h(x, Rfast::Dist::wave_hedges, parallel);
	}
	else if (method == "motyka")
	{
		return DistTotal::dist_h(x, Rfast::Dist::motyka, parallel);
	}
	else if (method == "harmonic_mean")
	{
		return DistTotal::dist_h(x, Rfast::Dist::harmonic_mean, parallel);
	}
	else if (method == "total_variation")
	{
		return DistTotal::dist_h(x, Rfast::Dist::total_variation, parallel);
	}
	else if (method == "sorensen")
	{
		return DistTotal::dist_h(x, Rfast::Dist::sorensen, parallel);
	}
	else if (method == "maximum")
	{
		return DistTotal::dist_h(x, Rfast::Dist::max, parallel);
	}
	else if (method == "minimum")
	{
		return DistTotal::dist_h(x, Rfast::Dist::min, parallel);
	}
	else if (method == "hellinger")
	{
		return sqr ? DistTotal::dist_h(x, 0.5, Rfast::Dist::hellinger<true>, parallel) : DistTotal::dist_h(x, 1.0 / std::sqrt(2.0), Rfast::Dist::hellinger<false>, parallel);
	}
	else if (method == "gower")
	{
		return DistTotal::dist_h(x, 1.0 / x.nrow(), Rfast::Dist::gower, parallel);
	}
	stop("Unsupported Method: %s", method);
}

RcppExport SEXP Rfast_total_dists(SEXP xSEXP, SEXP methodSEXP, SEXP sqrSEXP, SEXP pSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type sqr(sqrSEXP);
	traits::input_parameter<const int>::type p(pSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = total_dist(x, method, sqr, p, parallel);
	return __result;
	END_RCPP
}