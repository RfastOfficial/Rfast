// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "Rfast.h"
#include "mn.h"

using namespace arma;
using namespace Rcpp;
using std::string;

void euclidean_dista(mat &xnew, mat &x, mat &disa, const bool sqr)
{
	if (sqr)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = sum(square(x.each_col() - xnew.col(i)), 0).t();
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

void manhattan_dista(mat &xnew, mat &x, mat &disa)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = sum(square(x.each_col() - xnew.col(i)), 0).t() * 0.5;
	}
}

void hellinger_dista(mat &xnew, mat &x, mat &disa, const bool sqr)
{
	if (sqr)
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = sum(square(x.each_col() - xnew.col(i)), 0).t() * 0.5;
		}
	}
	else
	{
		constexpr double p = 1.0 / std::sqrt(2.0);
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			disa.col(i) = foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)).t() * p;
		}
	}
}

void max_dista(mat &xnew, mat &x, mat &disa)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = max(abs(x.each_col() - xnew.col(i)), 0).t();
	}
}

void min_dista(mat &xnew, mat &x, mat &disa)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = min(abs(x.each_col() - xnew.col(i)), 0).t();
	}
}

void minkowski_dista(mat &xnew, mat &x, mat &disa, const double p)
{
	const double p_1 = 1.0 / p;
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = pow(sum(pow(abs(x.each_col() - xnew.col(i)), p), 0), p_1).t();
	}
}

void canberra_dista(mat &xnew, mat &x, mat &disa)
{
	mat x_abs = abs(x);
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = sum(abs(x.each_col() - xnew.col(i)) / (x_abs.each_col() + abs(xnew.col(i))), 0).t();
	}
}

void total_variation_dista(mat &xnew, mat &x, mat &disa)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = sum(abs(x.each_col() - xnew.col(i)), 0).t() * 0.5;
	}
}

void kullback_leibler_dista(mat &xnew, mat &x, mat &disa, const bool parallel = false)
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
			disa.col(i) = colsum_with_condition<rowvec, std::isfinite>(m).t();
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
			disa.col(i) = colsum_with_condition<rowvec, std::isfinite>(m).t();
		}
	}
}

static bool check_if_is_finite(double x)
{
	return x > 0 and !R_IsNA(x);
}

void jensen_shannon_dista(mat &xnew, mat &x, mat &disa, const bool parallel = false)
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
			disa.col(i) = colsum_with_condition<rowvec, check_if_is_finite>(m).t();
		}
	}
	else
	{
		for (unsigned int i = 0; i < disa.n_cols; ++i)
		{
			mat xcolj = x.each_col() + xnew.col(i);
			mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
			mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
			disa.col(i) = colsum_with_condition<rowvec, check_if_is_finite>(m).t();
		}
	}
}

void bhattacharyya_dista(mat &xnew, mat &x, mat &disa)
{
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = sum(sqrt(x.each_col() % xnew.col(i)), 0).t();
	}
}

void itakura_saito_dista(mat &xnew, mat &x, mat &disa)
{
	mat log_x(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
	fill_with<std::log, double *, double *>(x.begin(), x.end(), log_x.begin());
	fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());
	for (unsigned int i = 0; i < disa.n_cols; ++i)
	{
		disa.col(i) = sum(x.each_col() / xnew.col(i) - (log_x.each_col() - log_xnew.col(i)) - 1, 0).t();
	}
}

//[[Rcpp::export]]
NumericMatrix dista(NumericMatrix Xnew, NumericMatrix X, const bool sqr = false, const string method = "", const double p = 0.0, const bool parallel = false)
{
	const int n = X.ncol(), nu = Xnew.ncol();
	mat xnew(Xnew.begin(), Xnew.nrow(), nu, false), x(X.begin(), X.nrow(), n, false);
	NumericMatrix disaa(n, nu);
	mat disa(disaa.begin(), n, nu, false);
	if (method == "euclidean")
	{
		euclidean_dista(xnew, x, disa, sqr);
	}
	else if (method == "manhattan")
	{
		manhattan_dista(xnew, x, disa);
	}
	else if (method == "hellinger")
	{
		hellinger_dista(xnew, x, disa, sqr);
	}
	else if (method == "maximum")
	{
		max_dista(xnew, x, disa);
	}
	else if (method == "minimum")
	{
		min_dista(xnew, x, disa);
	}
	else if (method == "minkowski")
	{
		minkowski_dista(xnew, x, disa, p);
	}
	else if (method == "canberra")
	{
		canberra_dista(xnew, x, disa);
	}
	else if (method == "bhattacharyya")
	{
		bhattacharyya_dista(xnew, x, disa);
	}
	else if (method == "jensen_shannon")
	{
		jensen_shannon_dista(xnew, x, disa, parallel);
	}
	else if (method == "itakura_saito")
	{
		itakura_saito_dista(xnew, x, disa);
	}
	else if (method == "total_variation")
	{
		total_variation_dista(xnew, x, disa);
	}
	else if (method == "kullback_leibler")
	{
		kullback_leibler_dista(xnew, x, disa, parallel);
	}
	else
		stop("Unsupported Method: %s", method);
	return disaa;
}

RcppExport SEXP Rfast_dista(SEXP XnewSEXP, SEXP XSEXP, SEXP sqrSEXP, SEXP methodSEXP, SEXP pSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type Xnew(XnewSEXP);
	traits::input_parameter<NumericMatrix>::type X(XSEXP);
	traits::input_parameter<const bool>::type sqr(sqrSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const int>::type p(pSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = dista(Xnew, X, sqr, method, p, parallel);
	return __result;

	END_RCPP
}

//[[Rcpp::export]]
IntegerMatrix dista_index(NumericMatrix Xnew, NumericMatrix X, const int k, const string type)
{
	const int n = X.ncol(), nu = Xnew.ncol(), p = Xnew.nrow();
	mat xnew(Xnew.begin(), p, nu, false), x(X.begin(), p, n, false);
	IntegerMatrix disaa(k, nu);
	imat disa(disaa.begin(), k, nu, false);
	if (type == "euclidean")
	{
		for (int i = 0; i < nu; ++i)
			disa.col(i) = get_k_indices(sum(square(x.each_col() - xnew.col(i)), 0), k);
	}
	else if (type == "manhattan")
	{
		for (int i = 0; i < nu; ++i)
			disa.col(i) = get_k_indices(sum(abs(x.each_col() - xnew.col(i)), 0), k);
	}
	else
		stop("Unknown type argument. you have to enter \"euclidean\" or \"manhattan\".");
	return disaa;
}

RcppExport SEXP Rfast_dista_index(SEXP XnewSEXP, SEXP XSEXP, SEXP kSEXP, SEXP methodSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type Xnew(XnewSEXP);
	traits::input_parameter<NumericMatrix>::type X(XSEXP);
	traits::input_parameter<const int>::type k(kSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	__result = dista_index(Xnew, X, k, method);
	return __result;

	END_RCPP
}

//[[Rcpp::export]]
colvec get_k_values(rowvec x, const int &k)
{
	sort(x.begin(), x.end());
	return conv_to<colvec>::from(x.subvec(0, k - 1));
}

//[[Rcpp::export]]
NumericMatrix dista_values(NumericMatrix Xnew, NumericMatrix X, const int k, const bool sqr, const string type)
{
	const int n = X.ncol(), nu = Xnew.ncol();
	mat xnew(Xnew.begin(), Xnew.nrow(), nu, false), x(X.begin(), X.nrow(), n, false);
	NumericMatrix disaa(k, nu);
	mat disa(disaa.begin(), k, nu, false);
	if (type == "euclidean")
	{
		if (sqr)
		{
			for (int i = 0; i < nu; ++i)
				disa.col(i) = get_k_values(sum(square(x.each_col() - xnew.col(i)), 0), k);
		}
		else
		{
			for (int i = 0; i < nu; ++i)
			{
				disa.col(i) = foreach<std::sqrt, colvec>(get_k_values(sum(square(x.each_col() - xnew.col(i)), 0), k));
			}
		}
	}
	else if (type == "manhattan")
	{
		for (int i = 0; i < nu; ++i)
			disa.col(i) = get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0), k);
	}
	else
		stop("Unknown type argument. you have to enter \"euclidean\" or \"manhattan\".");
	return disaa;
}

RcppExport SEXP Rfast_dista_values(SEXP XnewSEXP, SEXP XSEXP, SEXP kSEXP, SEXP sqrSEXP, SEXP methodSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type Xnew(XnewSEXP);
	traits::input_parameter<NumericMatrix>::type X(XSEXP);
	traits::input_parameter<const int>::type k(kSEXP);
	traits::input_parameter<const bool>::type sqr(sqrSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	__result = dista_values(Xnew, X, k, sqr, method);
	return __result;

	END_RCPP
}