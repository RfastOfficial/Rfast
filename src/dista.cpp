// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "Rfast.h"
#include "mn.h"

using namespace arma;
using namespace Rcpp;
using std::string;

namespace Dista
{

	void euclidean(mat &xnew, mat &x, mat &disa, const bool sqr, const unsigned int k)
	{
		if (sqr)
		{
			if (k > 0)
			{

				for (unsigned int i = 0; i < disa.n_cols; ++i)
				{
					disa.col(i) = get_k_values(sum(square(x.each_col() - xnew.col(i)), 0), k);
				}
			}
			else
			{
				for (unsigned int i = 0; i < disa.n_cols; ++i)
				{
					disa.col(i) = sum(square(x.each_col() - xnew.col(i)), 0).t();
				}
			}
		}
		else
		{
			if (k > 0)
			{

				for (unsigned int i = 0; i < disa.n_cols; ++i)
				{
					disa.col(i) = get_k_values(foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)), k);
				}
			}
			else
			{
				for (unsigned int i = 0; i < disa.n_cols; ++i)
				{
					disa.col(i) = foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)).t();
				}
			}
		}
	}

	void manhattan(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{

			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(abs(x.each_col() - xnew.col(i)), 0).t();
			}
		}
	}

	void sorensen(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(abs(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0).t();
			}
		}
	}

	void chi_square(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{

			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(square(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(square(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0).t();
			}
		}
	}

	void cosine(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		colvec norm_xnew = euclidean_norm(xnew).t();
		rowvec norm_x = euclidean_norm(x);
		if (k > 0)
		{

			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(x.each_col() % xnew.col(i), 0) / (norm_x * norm_xnew[i]), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = (sum(x.each_col() % xnew.col(i), 0) / (norm_x * norm_xnew[i])).t();
			}
		}
	}

	void soergel(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0) / colSumMaxs<rowvec>(x, xnew.col(i)), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(abs(x.each_col() - xnew.col(i)), 0).t() / colSumMaxs<colvec>(x, xnew.col(i));
			}
		}
	}

	void kulczynski(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0) / colSumMins<rowvec>(x, xnew.col(i)), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(abs(x.each_col() - xnew.col(i)), 0).t() / colSumMins<colvec>(x, xnew.col(i));
			}
		}
	}

	void motyka(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(1.0 - colSumMins<rowvec>(x, xnew.col(i)) / sum(abs(x.each_col() + xnew.col(i)), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = 1.0 - colSumMins<colvec>(x, xnew.col(i)) / sum(abs(x.each_col() + xnew.col(i)), 0).t();
			}
		}
	}

	void harmonic_mean(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(x.each_col() % xnew.col(i), 0) / sum(x.each_col() + xnew.col(i), 0), k) * 2.0;
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = (sum(x.each_col() % xnew.col(i), 0) / sum(x.each_col() + xnew.col(i), 0)).t() * 2.0;
			}
		}
	}

	void hellinger(mat &xnew, mat &x, mat &disa, const bool sqr, const unsigned int k)
	{
		if (sqr)
		{
			if (k > 0)
			{
				for (unsigned int i = 0; i < disa.n_cols; ++i)
				{
					disa.col(i) = get_k_values(sum(square(x.each_col() - xnew.col(i)), 0), k) * 0.5;
				}
			}
			else
			{
				for (unsigned int i = 0; i < disa.n_cols; ++i)
				{
					disa.col(i) = sum(square(x.each_col() - xnew.col(i)), 0).t() * 0.5;
				}
			}
		}
		else
		{
			const double p = 1.0 / std::sqrt(2.0);
			if (k > 0)
			{
				for (unsigned int i = 0; i < disa.n_cols; ++i)
				{
					disa.col(i) = get_k_values(foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)), k) * p;
				}
			}
			else
			{
				for (unsigned int i = 0; i < disa.n_cols; ++i)
				{
					disa.col(i) = foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)).t() * p;
				}
			}
		}
	}

	void max(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(max(abs(x.each_col() - xnew.col(i)), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = max(abs(x.each_col() - xnew.col(i)), 0).t();
			}
		}
	}

	void min(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{

			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(min(abs(x.each_col() - xnew.col(i)), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = min(abs(x.each_col() - xnew.col(i)), 0).t();
			}
		}
	}

	void minkowski(mat &xnew, mat &x, mat &disa, const double p, const unsigned int k)
	{
		const double p_1 = 1.0 / p;

		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(pow(sum(pow(abs(x.each_col() - xnew.col(i)), p), 0), p_1), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = pow(sum(pow(abs(x.each_col() - xnew.col(i)), p), 0), p_1).t();
			}
		}
	}

	void canberra(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		mat x_abs = abs(x);

		if (k > 0)
		{

			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)) / (x_abs.each_col() + abs(xnew.col(i))), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(abs(x.each_col() - xnew.col(i)) / (x_abs.each_col() + abs(xnew.col(i))), 0).t();
			}
		}
	}

	void total_variation(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{

			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0), k) * 0.5;
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(abs(x.each_col() - xnew.col(i)), 0).t() * 0.5;
			}
		}
	}

	void kullback_leibler(mat &xnew, mat &x, mat &disa, const unsigned int k, const bool parallel = false)
	{
		mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
		fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

		if (k > 0)
		{
			// #pragma omp parallel for if (parallel)
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
				disa.col(i) = get_k_values(colsum_with_condition<rowvec, std::isfinite>(m), k);
			}
		}
		else
		{
#pragma omp parallel for if (parallel)
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
				disa.col(i) = colsum_with_condition<colvec, std::isfinite>(m);
			}
		}
	}

	void jensen_shannon(mat &xnew, mat &x, mat &disa, const unsigned int k, const bool parallel = false)
	{
		mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
		const double log2 = std::log(2);
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
		fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());
		mat x_mod_log_xx = x % log_xx;

		if (k > 0)
		{
#pragma omp parallel for if (parallel)
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat xcolj = x.each_col() + xnew.col(i);
				mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
				mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
				disa.col(i) = get_k_values(colsum_with_condition<rowvec, check_if_is_finite>(m), k);
			}
		}
		else
		{
#pragma omp parallel for if (parallel)
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat xcolj = x.each_col() + xnew.col(i);
				mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
				mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
				disa.col(i) = colsum_with_condition<colvec, check_if_is_finite>(m);
			}
		}
	}

	void bhattacharyya(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{

			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(-log(sum(sqrt(x.each_col() % xnew.col(i)), 0)), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = -log(sum(sqrt(x.each_col() % xnew.col(i)), 0)).t();
			}
		}
	}

	void jeffries_matusita(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{

			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sqrt(2.0 - 2.0 * sum(sqrt(x.each_col() % xnew.col(i)), 0)), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sqrt(2.0 - 2.0 * sum(sqrt(x.each_col() % xnew.col(i)), 0)).t();
			}
		}
	}

	void itakura_saito(mat &xnew, mat &x, mat &disa, const unsigned int k, const bool parallel = false)
	{
		mat log_x(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_x.begin());
		fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

		if (k > 0)
		{
#pragma omp parallel for if (parallel)
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
			}
		}
		else
		{
#pragma omp parallel for if (parallel)
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
			}
		}
	}

	void wave_hedges(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)) / colMaxElems(x, xnew.col(i)), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(abs(x.each_col() - xnew.col(i)) / colMaxElems(x, xnew.col(i)), 0).t();
			}
		}
	}

	void gower(mat &xnew, mat &x, mat &disa, const unsigned int k)
	{
		const double p = 1.0 / x.n_rows;
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)) * p, 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = sum(abs(x.each_col() - xnew.col(i)) * p, 0).t();
			}
		}
	}

}

//[[Rcpp::export]]
NumericMatrix dista(NumericMatrix Xnew, NumericMatrix X, const string method = "", const bool sqr = false, const double p = 0.0, const unsigned int k = 0, const bool parallel = false)
{
	// if k is greater than 0 then rows are k size
	const int n = k > 0 ? k : X.ncol(), nu = Xnew.ncol();
	mat xnew(Xnew.begin(), Xnew.nrow(), nu, false), x(X.begin(), X.nrow(), X.ncol(), false);
	NumericMatrix disaa(n, nu);
	mat disa(disaa.begin(), n, nu, false);
	if (method == "euclidean")
	{
		Dista::euclidean(xnew, x, disa, sqr, k);
	}
	else if (method == "manhattan")
	{
		Dista::manhattan(xnew, x, disa, k);
	}
	else if (method == "hellinger")
	{
		Dista::hellinger(xnew, x, disa, sqr, k);
	}
	else if (method == "maximum")
	{
		Dista::max(xnew, x, disa, k);
	}
	else if (method == "minimum")
	{
		Dista::min(xnew, x, disa, k);
	}
	else if (method == "minkowski")
	{
		Dista::minkowski(xnew, x, disa, p, k);
	}
	else if (method == "canberra")
	{
		Dista::canberra(xnew, x, disa, k);
	}
	else if (method == "bhattacharyya")
	{
		Dista::bhattacharyya(xnew, x, disa, k);
	}
	else if (method == "jensen_shannon")
	{
		Dista::jensen_shannon(xnew, x, disa, k, parallel);
	}
	else if (method == "itakura_saito")
	{
		Dista::itakura_saito(xnew, x, disa, k, parallel);
	}
	else if (method == "total_variation")
	{
		Dista::total_variation(xnew, x, disa, k);
	}
	else if (method == "kullback_leibler")
	{
		Dista::kullback_leibler(xnew, x, disa, k, parallel);
	}
	else if (method == "chi_square")
	{
		Dista::chi_square(xnew, x, disa, k);
	}
	else if (method == "sorensen")
	{
		Dista::sorensen(xnew, x, disa, k);
	}
	else if (method == "soergel")
	{
		Dista::soergel(xnew, x, disa, k);
	}
	else if (method == "cosine")
	{
		Dista::cosine(xnew, x, disa, k);
	}
	else if (method == "wave_hedges")
	{
		Dista::wave_hedges(xnew, x, disa, k);
	}
	else if (method == "motyka")
	{
		Dista::motyka(xnew, x, disa, k);
	}
	else if (method == "harmonic_mean")
	{
		Dista::harmonic_mean(xnew, x, disa, k);
	}
	else if (method == "jeffries_matusita")
	{
		Dista::jeffries_matusita(xnew, x, disa, k);
	}
	else if (method == "gower")
	{
		Dista::gower(xnew, x, disa, k);
	}
	else if (method == "kulczynski")
	{
		Dista::kulczynski(xnew, x, disa, k);
	}
	else
		stop("Unsupported Method: %s", method);
	return disaa;
}

namespace DistaIndices
{

	void euclidean(mat &xnew, mat &x, imat &disa, const bool sqr, const unsigned int k)
	{
		if (sqr)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_indices(sum(square(x.each_col() - xnew.col(i)), 0), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_indices(foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)), k);
			}
		}
	}

	void manhattan(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0), k);
		}
	}

	void sorensen(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0), k);
		}
	}

	void chi_square(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(square(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0), k);
		}
	}

	void hellinger(mat &xnew, mat &x, imat &disa, const bool sqr, const unsigned int k)
	{
		if (sqr)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_indices(sum(square(x.each_col() - xnew.col(i)), 0) * 0.5, k);
			}
		}
		else
		{
			const double p = 1.0 / std::sqrt(2.0);
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_indices(foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)) * p, k);
			}
		}
	}

	void max(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(max(abs(x.each_col() - xnew.col(i)), 0), k);
		}
	}

	void min(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(min(abs(x.each_col() - xnew.col(i)), 0), k);
		}
	}

	void minkowski(mat &xnew, mat &x, imat &disa, const double p, const unsigned int k)
	{
		const double p_1 = 1.0 / p;

		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(pow(sum(pow(abs(x.each_col() - xnew.col(i)), p), 0), p_1), k);
		}
	}

	void canberra(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		mat x_abs = abs(x);

		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)) / (x_abs.each_col() + abs(xnew.col(i))), 0), k);
		}
	}

	void total_variation(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0) * 0.5, k);
		}
	}

	void soergel(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0) / colSumMaxs<colvec>(x, xnew.col(i)), k);
		}
	}

	void kulczynski(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0) / colSumMins<colvec>(x, xnew.col(i)), k);
		}
	}

	void kullback_leibler(mat &xnew, mat &x, imat &disa, const unsigned int k, const bool parallel = false)
	{
		mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
		fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

#pragma omp parallel for if (parallel)
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, std::isfinite>(m), k);
		}
	}

	void jensen_shannon(mat &xnew, mat &x, imat &disa, const unsigned int k, const bool parallel = false)
	{
		mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
		const double log2 = std::log(2);
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
		fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());
		mat x_mod_log_xx = x % log_xx;

#pragma omp parallel for if (parallel)
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat xcolj = x.each_col() + xnew.col(i);
			mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
			mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, check_if_is_finite>(m), k);
		}
	}

	void bhattacharyya(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(-log(sum(sqrt(x.each_col() % xnew.col(i)), 0)), k);
		}
	}

	void cosine(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		colvec norm_xnew = euclidean_norm(xnew).t();
		rowvec norm_x = euclidean_norm(x);
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(x.each_col() % xnew.col(i), 0) / (norm_x * norm_xnew[i]), k);
		}
	}

	void wave_hedges(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)) / colMaxElems(x, xnew.col(i)), 0), k);
		}
	}

	void motyka(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(1.0 - colSumMins<rowvec>(x, xnew.col(i)) / sum(abs(x.each_col() + xnew.col(i)), 0), k);
		}
	}

	void harmonic_mean(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(x.each_col() % xnew.col(i), 0) / sum(x.each_col() + xnew.col(i), 0) * 2.0, k);
		}
	}

	void jeffries_matusita(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sqrt(2.0 - 2.0 * sum(sqrt(x.each_col() % xnew.col(i)), 0)), k);
		}
	}

	void itakura_saito(mat &xnew, mat &x, imat &disa, const unsigned int k, const bool parallel = false)
	{
		mat log_x(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
		fill_with<std::log, double *, double *>(x.begin(), x.end(), log_x.begin());
		fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

#pragma omp parallel for if (parallel)
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat m = (x.each_col() / xnew.col(i)) - (log_x.each_col() - log_xnew.col(i)) - 1;
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, std::isfinite>(m), k);
		}
	}

	void gower(mat &xnew, mat &x, imat &disa, const unsigned int k)
	{
		const double p = 1.0 / x.n_rows;
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)) * p, 0), k);
		}
	}
}

IntegerMatrix dista_index(NumericMatrix Xnew, NumericMatrix X, const string method = "", const bool sqr = false, const double p = 0.0, const unsigned int k = 0, const bool parallel = false)
{
	// if k is greater than 0 then rows are k size
	const int n = k > 0 ? k : X.ncol(), nu = Xnew.ncol();
	mat xnew(Xnew.begin(), Xnew.nrow(), nu, false), x(X.begin(), X.nrow(), X.ncol(), false);
	IntegerMatrix disaa(n, nu);
	imat disa(disaa.begin(), n, nu, false);
	if (method == "euclidean")
	{
		DistaIndices::euclidean(xnew, x, disa, sqr, k);
	}
	else if (method == "manhattan")
	{
		DistaIndices::manhattan(xnew, x, disa, k);
	}
	else if (method == "hellinger")
	{
		DistaIndices::hellinger(xnew, x, disa, sqr, k);
	}
	else if (method == "maximum")
	{
		DistaIndices::max(xnew, x, disa, k);
	}
	else if (method == "minimum")
	{
		DistaIndices::min(xnew, x, disa, k);
	}
	else if (method == "minkowski")
	{
		DistaIndices::minkowski(xnew, x, disa, p, k);
	}
	else if (method == "canberra")
	{
		DistaIndices::canberra(xnew, x, disa, k);
	}
	else if (method == "bhattacharyya")
	{
		DistaIndices::bhattacharyya(xnew, x, disa, k);
	}
	else if (method == "jensen_shannon")
	{
		DistaIndices::jensen_shannon(xnew, x, disa, k, parallel);
	}
	else if (method == "itakura_saito")
	{
		DistaIndices::itakura_saito(xnew, x, disa, k, parallel);
	}
	else if (method == "total_variation")
	{
		DistaIndices::total_variation(xnew, x, disa, k);
	}
	else if (method == "kullback_leibler")
	{
		DistaIndices::kullback_leibler(xnew, x, disa, k, parallel);
	}
	else if (method == "chi_square")
	{
		DistaIndices::chi_square(xnew, x, disa, k);
	}
	else if (method == "sorensen")
	{
		DistaIndices::sorensen(xnew, x, disa, k);
	}
	else if (method == "soergel")
	{
		DistaIndices::soergel(xnew, x, disa, k);
	}
	else if (method == "cosine")
	{
		DistaIndices::cosine(xnew, x, disa, k);
	}
	else if (method == "wave_hedges")
	{
		DistaIndices::wave_hedges(xnew, x, disa, k);
	}
	else if (method == "motyka")
	{
		DistaIndices::motyka(xnew, x, disa, k);
	}
	else if (method == "harmonic_mean")
	{
		DistaIndices::harmonic_mean(xnew, x, disa, k);
	}
	else if (method == "jeffries_matusita")
	{
		DistaIndices::jeffries_matusita(xnew, x, disa, k);
	}
	else if (method == "gower")
	{
		DistaIndices::gower(xnew, x, disa, k);
	}
	else if (method == "kulczynski")
	{
		DistaIndices::kulczynski(xnew, x, disa, k);
	}
	else
		stop("Unsupported Method: %s", method);
	return disaa;
}

RcppExport SEXP Rfast_dista(SEXP XnewSEXP, SEXP XSEXP, SEXP methodSEXP, SEXP sqrSEXP, SEXP pSEXP, SEXP kSEXP, SEXP indexSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type Xnew(XnewSEXP);
	traits::input_parameter<NumericMatrix>::type X(XSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type sqr(sqrSEXP);
	traits::input_parameter<const double>::type p(pSEXP);
	traits::input_parameter<const unsigned int>::type k(kSEXP);
	traits::input_parameter<const bool>::type index(indexSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);

	if (index)
		__result = dista_index(Xnew, X, method, sqr, p, k, parallel);
	else
		__result = dista(Xnew, X, method, sqr, p, k, parallel);
	return __result;

	END_RCPP
}
