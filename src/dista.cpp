// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "Rfast.h"
#include "mn.h"

using namespace arma;
using namespace Rcpp;
using std::string;

void euclidean_dista(mat &xnew, mat &x, mat &disa, const bool sqr, const unsigned int k)
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

void manhattan_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
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

void sorensen_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
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

void chi_square_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
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

void cosine_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
{
	colvec norm_xnew = euclidean_norm(xnew), norm_x = euclidean_norm(x);
	if (k > 0)
	{

		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_values(sum(x.each_col() % xnew.col(i), 0).t() / (norm_x * norm_xnew[i]), k);
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = sum(x.each_col() % xnew.col(i), 0).t() / (norm_x * norm_xnew[i]);
		}
	}
}

void soergel_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
{
	if (k > 0)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0) / colSumMaxs<rowvec>(x,xnew.col(i)), k);
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = sum(abs(x.each_col() - xnew.col(i)), 0).t() / colSumMaxs<colvec>(x,xnew.col(i));
		}
	}
}

void hellinger_dista(mat &xnew, mat &x, mat &disa, const bool sqr, const unsigned int k)
{
	if (sqr)
	{
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(sum(square(x.each_col() - xnew.col(i)), 0) * 0.5, k);
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
		constexpr double p = 1.0 / std::sqrt(2.0);
		if (k > 0)
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				disa.col(i) = get_k_values(foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)) * p, k);
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

void max_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
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

void min_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
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

void minkowski_dista(mat &xnew, mat &x, mat &disa, const double p, const unsigned int k)
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

void canberra_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
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

void total_variation_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
{
	if (k > 0)
	{

		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0) * 0.5, k);
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

void kullback_leibler_dista(mat &xnew, mat &x, mat &disa, const unsigned int k, const bool parallel = false)
{
	mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
	fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
	fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

	if (k > 0)
	{
		if (parallel)
		{
#pragma omp parallel for
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
				disa.col(i) = get_k_values(colsum_with_condition<colvec, std::isfinite>(m), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
				disa.col(i) = get_k_values(colsum_with_condition<colvec, std::isfinite>(m), k);
			}
		}
	}
	else
	{
		if (parallel)
		{
#pragma omp parallel for
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
				disa.col(i) = colsum_with_condition<colvec, std::isfinite>(m);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
				disa.col(i) = colsum_with_condition<colvec, std::isfinite>(m);
			}
		}
	}
}

void jensen_shannon_dista(mat &xnew, mat &x, mat &disa, const unsigned int k, const bool parallel = false)
{
	mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
	constexpr double log2 = std::log(2);
	fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
	fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());
	mat x_mod_log_xx = x % log_xx;

	if (k > 0)
	{
		if (parallel)
		{
#pragma omp parallel for
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat xcolj = x.each_col() + xnew.col(i);
				mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
				mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
				disa.col(i) = get_k_values(colsum_with_condition<colvec, check_if_is_finite>(m), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat xcolj = x.each_col() + xnew.col(i);
				mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
				mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
				disa.col(i) = get_k_values(colsum_with_condition<colvec, check_if_is_finite>(m), k);
			}
		}
	}
	else
	{
		if (parallel)
		{
#pragma omp parallel for
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat xcolj = x.each_col() + xnew.col(i);
				mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
				mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
				disa.col(i) = colsum_with_condition<colvec, check_if_is_finite>(m);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat xcolj = x.each_col() + xnew.col(i);
				mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
				mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
				disa.col(i) = colsum_with_condition<colvec, check_if_is_finite>(m);
			}
		}
	}
}

void bhattacharyya_dista(mat &xnew, mat &x, mat &disa, const unsigned int k)
{
	if (k > 0)
	{

		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_values(sum(sqrt(x.each_col() % xnew.col(i)), 0), k);
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = sum(sqrt(x.each_col() % xnew.col(i)), 0).t();
		}
	}
}

void itakura_saito_dista(mat &xnew, mat &x, mat &disa, const unsigned int k, const bool parallel = false)
{
	mat log_x(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
	fill_with<std::log, double *, double *>(x.begin(), x.end(), log_x.begin());
	fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

	if (k > 0)
	{
		if (parallel)
		{
#pragma omp parallel for
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_x.each_col() - log_xnew.col(i));
				disa.col(i) = get_k_values(colsum_with_condition<colvec, std::isfinite>(m), k);
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_x.each_col() - log_xnew.col(i));
				disa.col(i) = get_k_values(colsum_with_condition<colvec, std::isfinite>(m), k);
			}
		}
	}
	else
	{
		if (parallel)
		{
#pragma omp parallel for
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_x.each_col() - log_xnew.col(i));
				disa.col(i) = colsum_with_condition<colvec, std::isfinite>(m).t();
			}
		}
		else
		{
			for (unsigned int i = 0; i < disa.n_cols; ++i)
			{
				mat m = (x.each_col() - xnew.col(i)) % (log_x.each_col() - log_xnew.col(i));
				disa.col(i) = colsum_with_condition<colvec, std::isfinite>(m).t();
			}
		}
	}
}

void wave_hedges(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
  	mat x_max = max(x,0),xnew_max = max(xnew,0);
	if (k > 0)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0).t() / max(x_max[i],xnew_max[i]), k);
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = sum(abs(x.each_col() - xnew.col(i)), 0).t() / max(x_max[i],xnew_max[i]).t();
		}
	}
}

//[[Rcpp::export]]
NumericMatrix dista(NumericMatrix Xnew, NumericMatrix X, const string method = "", const bool sqr = false, const double p = 0.0, const unsigned int k = 0, const bool parallel = false)
{
	const int n = X.ncol(), nu = Xnew.ncol();
	mat xnew(Xnew.begin(), Xnew.nrow(), nu, false), x(X.begin(), X.nrow(), n, false);
	NumericMatrix disaa(n, nu);
	mat disa(disaa.begin(), n, nu, false);
	if (method == "euclidean")
	{
		euclidean_dista(xnew, x, disa, sqr, k);
	}
	else if (method == "manhattan")
	{
		manhattan_dista(xnew, x, disa, k);
	}
	else if (method == "hellinger")
	{
		hellinger_dista(xnew, x, disa, sqr, k);
	}
	else if (method == "maximum")
	{
		max_dista(xnew, x, disa, k);
	}
	else if (method == "minimum")
	{
		min_dista(xnew, x, disa, k);
	}
	else if (method == "minkowski")
	{
		minkowski_dista(xnew, x, disa, p, k);
	}
	else if (method == "canberra")
	{
		canberra_dista(xnew, x, disa, k);
	}
	else if (method == "bhattacharyya")
	{
		bhattacharyya_dista(xnew, x, disa, k);
	}
	else if (method == "jensen_shannon")
	{
		jensen_shannon_dista(xnew, x, disa, k, parallel);
	}
	else if (method == "itakura_saito")
	{
		itakura_saito_dista(xnew, x, disa, k, parallel);
	}
	else if (method == "total_variation")
	{
		total_variation_dista(xnew, x, disa, k);
	}
	else if (method == "kullback_leibler")
	{
		kullback_leibler_dista(xnew, x, disa, k, parallel);
	}
	else if (method == "chi_square")
	{
		chi_square_dista(xnew, x, disa, k);
	}
	else if (method == "sorensen")
	{
		sorensen_dista(xnew, x, disa, k);
	}
	else if (method == "soergel")
	{
		soergel_dista(xnew, x, disa, k);
	}
	else if (method == "cosine")
	{
		cosine_dista(xnew, x, disa, k);
	}
	else if (method == "wave_hedges")
	{
		return wave_hedges_dista(x);
	}
	else
		stop("Unsupported Method: %s", method);
	return disaa;
}

void euclidean_dista_indices(mat &xnew, mat &x, imat &disa, const bool sqr, const unsigned int k)
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

void manhattan_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0), k);
	}
}

void sorensen_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0), k);
	}
}

void chi_square_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(square(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0), k);
	}
}

void hellinger_dista_indices(mat &xnew, mat &x, imat &disa, const bool sqr, const unsigned int k)
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
		constexpr double p = 1.0 / std::sqrt(2.0);
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = get_k_indices(foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)) * p, k);
		}
	}
}

void max_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(max(abs(x.each_col() - xnew.col(i)), 0), k);
	}
}

void min_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(min(abs(x.each_col() - xnew.col(i)), 0), k);
	}
}

void minkowski_dista_indices(mat &xnew, mat &x, imat &disa, const double p, const unsigned int k)
{
	const double p_1 = 1.0 / p;

	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(pow(sum(pow(abs(x.each_col() - xnew.col(i)), p), 0), p_1), k);
	}
}

void canberra_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	mat x_abs = abs(x);

	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)) / (x_abs.each_col() + abs(xnew.col(i))), 0), k);
	}
}

void total_variation_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0) * 0.5, k);
	}
}

void soergel_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0) / colSumMaxs<colvec>(x,xnew.col(i)), k);
	}
}

void kullback_leibler_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k, const bool parallel = false)
{
	mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
	fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
	fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

	if (parallel)
	{
#pragma omp parallel for
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, std::isfinite>(m), k);
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, std::isfinite>(m), k);
		}
	}
}

void jensen_shannon_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k, const bool parallel = false)
{
	mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
	constexpr double log2 = std::log(2);
	fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
	fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());
	mat x_mod_log_xx = x % log_xx;

	if (parallel)
	{
#pragma omp parallel for
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat xcolj = x.each_col() + xnew.col(i);
			mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
			mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, check_if_is_finite>(m), k);
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat xcolj = x.each_col() + xnew.col(i);
			mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
			mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, check_if_is_finite>(m), k);
		}
	}
}

void bhattacharyya_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(sqrt(x.each_col() % xnew.col(i)), 0), k);
	}
}

void cosine_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
	colvec norm_xnew = euclidean_norm(xnew), norm_x = euclidean_norm(x);
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(x.each_col() % xnew.col(i), 0).t() / (norm_x * norm_xnew[i]), k);
	}
}

void wave_hedges_indices(mat &xnew, mat &x, imat &disa, const unsigned int k)
{
  	mat x_max = max(x,0),xnew_max = max(xnew,0);
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0).t() / max(x_max[i],xnew_max[i]), k);
	}
}

void itakura_saito_dista_indices(mat &xnew, mat &x, imat &disa, const unsigned int k, const bool parallel = false)
{
	mat log_x(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
	fill_with<std::log, double *, double *>(x.begin(), x.end(), log_x.begin());
	fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

	if (parallel)
	{
#pragma omp parallel for
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat m = (x.each_col() - xnew.col(i)) % (log_x.each_col() - log_xnew.col(i));
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, std::isfinite>(m), k);
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat m = (x.each_col() - xnew.col(i)) % (log_x.each_col() - log_xnew.col(i));
			disa.col(i) = get_k_indices(colsum_with_condition<colvec, std::isfinite>(m), k);
		}
	}
}

IntegerMatrix dista_index(NumericMatrix Xnew, NumericMatrix X, const string method = "", const bool sqr = false, const double p = 0.0, const unsigned int k = 0, const bool parallel = false)
{
	const int n = X.ncol(), nu = Xnew.ncol();
	mat xnew(Xnew.begin(), Xnew.nrow(), nu, false), x(X.begin(), X.nrow(), n, false);
	IntegerMatrix disaa(n, nu);
	imat disa(disaa.begin(), n, nu, false);
	if (method == "euclidean")
	{
		euclidean_dista_indices(xnew, x, disa, sqr, k);
	}
	else if (method == "manhattan")
	{
		manhattan_dista_indices(xnew, x, disa, k);
	}
	else if (method == "hellinger")
	{
		hellinger_dista_indices(xnew, x, disa, sqr, k);
	}
	else if (method == "maximum")
	{
		max_dista_indices(xnew, x, disa, k);
	}
	else if (method == "minimum")
	{
		min_dista_indices(xnew, x, disa, k);
	}
	else if (method == "minkowski")
	{
		minkowski_dista_indices(xnew, x, disa, p, k);
	}
	else if (method == "canberra")
	{
		canberra_dista_indices(xnew, x, disa, k);
	}
	else if (method == "bhattacharyya")
	{
		bhattacharyya_dista_indices(xnew, x, disa, k);
	}
	else if (method == "jensen_shannon")
	{
		jensen_shannon_dista_indices(xnew, x, disa, k, parallel);
	}
	else if (method == "itakura_saito")
	{
		itakura_saito_dista_indices(xnew, x, disa, k, parallel);
	}
	else if (method == "total_variation")
	{
		total_variation_dista_indices(xnew, x, disa, k);
	}
	else if (method == "kullback_leibler")
	{
		kullback_leibler_dista_indices(xnew, x, disa, k, parallel);
	}
	else if (method == "chi_square")
	{
		chi_square_dista_indices(xnew, x, disa, k);
	}
	else if (method == "sorensen")
	{
		sorensen_dista_indices(xnew, x, disa, k);
	}
	else if (method == "soergel")
	{
		soergel_dista_indices(xnew, x, disa, k);
	}
	else if (method == "cosine")
	{
		cosine_dista_indices(xnew, x, disa, k);
	}
	else if (method == "wave_hedges")
	{
		return wave_hedges_dista_indices(x);
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
