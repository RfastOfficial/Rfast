
// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"
#include <chrono>
#include <random>
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
SEXP col_all(SEXP x)
{
	const int n = Rf_ncols(x), p = Rf_nrows(x);
	SEXP f = Rf_allocVector(LGLSXP, n);
	int *start = LOGICAL(x), *end = start + p, *ff = LOGICAL(f);
	for (int i = 0; i < n; ++i, ++ff)
	{
		*ff = my_all(start, end);
		start = end;
		end += p;
	}
	return f;
}

RcppExport SEXP Rfast_col_all(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_all(x);
	return __result;
	END_RCPP
}

LogicalVector row_all(LogicalMatrix x)
{
	const int n = x.nrow();
	LogicalVector f(n);
	for (int i = 0; i < n; ++i)
		f[i] = as<bool>(all(x.row(i)));
	return f;
}

RcppExport SEXP Rfast_row_all(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<LogicalMatrix>::type x(xSEXP);
	__result = row_all(x);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////

SEXP col_any(SEXP x)
{
	const int n = Rf_ncols(x), p = Rf_nrows(x);
	SEXP f = Rf_allocVector(LGLSXP, n);
	int *start = LOGICAL(x), *end = start + p, *ff = LOGICAL(f);
	for (int i = 0; i < n; ++i, ++ff)
	{
		*ff = my_any(start, end);
		start = end;
		end += p;
	}
	return f;
}

SEXP row_any(SEXP x)
{
	int nrow = Rf_nrows(x);
	SEXP F = PROTECT(Rf_allocVector(LGLSXP, nrow));
	int *xx = INTEGER(x), *endx = xx + LENGTH(x), *f = INTEGER(F), *startx, *startf;
	const int *endf = f + LENGTH(F);
	for (startf = f; startf != endf; ++startf)
		*startf = 0;
	while (xx != endx)
	{
		for (startf = f, startx = xx, xx += nrow; startx != xx; ++startf, ++startx)
		{
			if (*startx)
			{
				*startf = 1;
			}
		}
	}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_row_any(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = row_any(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_any(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_any(x);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////

IntegerVector col_count_values(NumericMatrix x, NumericVector values)
{
	const int n = values.size();
	IntegerVector f(n);
	for (int i = 0; i < n; ++i)
	{
		f[i] = count_value_helper<NumericVector, double>(x.column(i), values[i]);
	}
	return f;
}

RcppExport SEXP Rfast_col_count_values(SEXP xSEXP, SEXP valuesSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericVector>::type values(valuesSEXP);
	__result = col_count_values(x, values);
	return __result;
	END_RCPP
}

IntegerVector row_count_values(NumericMatrix x, NumericVector values)
{
	const int n = values.size();
	IntegerVector f(n);
	for (int i = 0; i < n; ++i)
	{
		f[i] = count_value_helper<NumericVector, double>(x.row(i), values[i]);
	}
	return f;
}

RcppExport SEXP Rfast_row_count_values(SEXP xSEXP, SEXP valuesSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericVector>::type values(valuesSEXP);
	__result = row_count_values(x, values);
	return __result;
	END_RCPP
}

//////////////////////////////////////////////////////

SEXP col_false(SEXP x)
{
	const int p = Rf_ncols(x), n = Rf_nrows(x);
	SEXP f = Rf_allocVector(INTSXP, p);
	int *ff = INTEGER(f), *xx = LOGICAL(x), *endx = xx + LENGTH(x);
	for (; xx != endx; xx += p, ++ff)
		*ff = n - True(xx, xx + p);
	return f;
}

SEXP row_false(SEXP x)
{
	int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	SEXP F = PROTECT(Rf_allocVector(INTSXP, nrow));
	int *xx = INTEGER(x), *end = xx + ncol * nrow, *f = INTEGER(F), *startx, *startf;
	const int *endf = f + LENGTH(F);
	for (startf = f; startf != endf; ++startf)
		*startf = ncol;
	while (xx != end)
	{
		for (startf = f, startx = xx, xx += nrow; startx != xx; ++startf, ++startx)
		{
			*startf -= *startx;
		}
	}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_row_false(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = row_false(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_false(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_false(x);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////

IntegerVector col_len_sort_un_int(IntegerMatrix x)
{
	const int p = x.ncol();
	IntegerVector f(p);
	for (int i = 0; i < p; ++i)
		f[i] = len_sort_unique_int(x.column(i));
	return f;
}

IntegerVector row_len_sort_un_int(IntegerMatrix x)
{
	const unsigned int p = x.nrow();
	IntegerVector F(p);
	IntegerVector::iterator FF = F.begin();
	for (int i = 0; FF != F.end(); ++FF, ++i)
	{
		*FF = len_sort_unique_int(x.row(i));
	}
	return F;
}

RcppExport SEXP Rfast_row_len_sort_un_int(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<IntegerMatrix>::type x(xSEXP);
	__result = row_len_sort_un_int(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_len_sort_un_int(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<IntegerMatrix>::type x(xSEXP);
	__result = col_len_sort_un_int(x);
	return __result;
	END_RCPP
}

//////////////////////////////////////////////////////////

RcppExport SEXP Rfast_col_mads(SEXP xSEXP, SEXP methodSEXP, SEXP na_rmSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = Rf_isMatrix(xSEXP) ? Rfast::colMads(NumericMatrix(xSEXP), method, na_rm, parallel) : Rfast::colMads(DataFrame(xSEXP), method, na_rm, parallel);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_row_mads(SEXP xSEXP, SEXP methodSEXP, SEXP na_rmSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = Rfast::rowMads(x, method, na_rm, parallel);
	return __result;
	END_RCPP
}

////////////////////////////////////////////////////

SEXP col_max_indices(NumericMatrix x)
{
	unsigned int i = 0, p = x.ncol();
	arma::mat X = arma::mat(x.begin(), x.nrow(), p, false);
	SEXP F = PROTECT(Rf_allocVector(INTSXP, p));
	int *FF = INTEGER(F);
	for (; i < p; ++i, ++FF)
		*FF = (X.col(i)).index_max() + 1;
	UNPROTECT(1);
	return F;
}

SEXP col_max(SEXP x, const bool parallel)
{
	int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	if (parallel)
	{
		NumericMatrix X(x);
		mat xx(X.begin(), nrow, ncol, false);
		NumericVector f(ncol);
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < ncol; ++i)
		{
			f[i] = *max_element(xx.begin_col(i), xx.end_col(i));
		}
		return f;
	}
	else
	{
		SEXP F;
		switch (Rfast::Type::type(x))
		{
		case Rfast::Type::Types::REAL:
		{
			F = PROTECT(Rf_allocVector(REALSXP, ncol));
			double *xx = REAL(x), *end = xx + ncol * nrow, *f = REAL(F);
			for (; xx != end; xx += nrow, ++f)
				maximum<double>(xx, xx + nrow, *f);
			break;
		}
		default:
		{
			F = PROTECT(Rf_allocVector(INTSXP, ncol));
			int *xx = INTEGER(x), *end = xx + ncol * nrow, *f = INTEGER(F);
			for (; xx != end; xx += nrow, ++f)
				maximum<int>(xx, xx + nrow, *f);
			break;
		}
		}
		UNPROTECT(1);
		return F;
	}
}

// find the maximum index of its collumn
RcppExport SEXP Rfast_col_max_indices(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = col_max_indices(x);
	return __result;
	END_RCPP
}

// find the maximum value of its collumn
RcppExport SEXP Rfast_col_max(SEXP x, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	if (Rf_isMatrix(x))
	{
		__result = col_max(x, parallel);
	}
	else
	{
		__result = Rfast::colMaxs(DataFrame(x), parallel);
	}
	return __result;
	END_RCPP
}

// [[Rcpp::export]]
SEXP row_max_indices(NumericMatrix x)
{
	const int p = x.nrow();
	mat X = mat(x.begin(), p, x.ncol(), false);
	SEXP F = PROTECT(Rf_allocVector(INTSXP, p));
	int *FF = INTEGER(F);
	for (int i = 0; i < p; ++i)
	{
		FF[i] = (X.row(i)).index_max() + 1;
	}
	UNPROTECT_PTR(F);
	return F;
}

SEXP row_max(SEXP x)
{
	int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	SEXP F;
	F = PROTECT(Rf_allocVector(REALSXP, nrow));
	double *xx = REAL(x), *end = xx + ncol * nrow, *f = REAL(F), *x3, *ff;
	const double *endf = f + LENGTH(F);
	for (ff = f; ff != endf; ++ff, ++xx)
		*ff = *xx;
	for (; xx != end;)
		for (ff = f, x3 = xx, xx += nrow; x3 != xx; ++ff, ++x3)
		{
			*ff = std::max(*ff, *x3);
		}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_row_max_indices(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = row_max_indices(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_row_max(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = row_max(x);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////////////

NumericVector col_means(DataFrame x, const bool parallel = false)
{
	int i = 0;
	NumericVector f(x.size());
	if (parallel)
	{
		colvec ff(f.begin(), f.size(), false);
#pragma omp parallel for
		for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
		{
			colvec y;
			int i;
#pragma omp critical
			{
				NumericVector yy;
				yy = *s;
				y = colvec(yy.begin(), yy.size(), false);
				i = s - x.begin();
			}
			ff[i] = mean(y);
		}
	}
	else
	{
		NumericVector y(x.nrows());
		for (auto c : x)
		{
			y = c;
			colvec yy(y.begin(), y.size(), false);
			f[i++] = mean(yy);
		}
	}
	f.names() = x.names();
	return f;
}

NumericVector col_means(NumericMatrix x, const bool parallel = false)
{
	mat xx;
	const int n = x.ncol();
	NumericVector F(n);
	rowvec ff(F.begin(), n, false);
	if (parallel)
	{
		xx = mat(x.begin(), x.nrow(), n, false);
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			ff[i] = mean(xx.col(i));
		}
		UNPROTECT(1);
	}
	else
	{
		xx = mat(x.begin(), x.nrow(), x.ncol(), false);
		ff = mean(xx, 0);
	}
	return F;
}

colvec row_means(NumericMatrix x)
{
	mat X = mat(x.begin(), x.nrow(), x.ncol(), false);
	return mean(X, 1);
}

RcppExport SEXP Rfast_row_means(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = row_means(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_means(SEXP xSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = Rf_isMatrix(xSEXP) ? col_means(NumericMatrix(xSEXP), parallel) : col_means(DataFrame(xSEXP), parallel);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////

// colMedians
RcppExport SEXP Rfast_col_meds(SEXP xSEXP, SEXP na_rmSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	if (Rf_isMatrix(xSEXP))
	{
		NumericMatrix x(xSEXP);
		__result = Rfast::colMedian(x, na_rm, parallel);
	}
	else
	{
		DataFrame x(xSEXP);
		__result = Rfast::colMedian(x, na_rm, parallel);
	}
	return __result;
	END_RCPP
}

// rowMedians
RcppExport SEXP Rfast_row_meds(SEXP xSEXP, SEXP na_rmSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = Rfast::rowMedian(x, na_rm, parallel);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////

SEXP col_min_indices(NumericMatrix x)
{
	unsigned int i = 0, p = x.ncol();
	mat X = mat(x.begin(), x.nrow(), p, false);
	SEXP F = PROTECT(Rf_allocVector(INTSXP, p));
	int *FF = INTEGER(F);
	for (; i < p; ++i, ++FF)
		*FF = (X.col(i)).index_min() + 1;
	UNPROTECT(1);
	return F;
}

SEXP col_min(SEXP x, const bool parallel)
{
	int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	if (parallel)
	{
		NumericMatrix X(x);
		mat xx(X.begin(), nrow, ncol, false);
		NumericVector f(ncol);
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < ncol; ++i)
		{
			f[i] = *min_element(xx.begin_col(i), xx.end_col(i));
		}
		return f;
	}
	else
	{
		SEXP F;
		switch (Rfast::Type::type(x))
		{
		case Rfast::Type::Types::REAL:
		{
			F = PROTECT(Rf_allocVector(REALSXP, ncol));
			double *xx = REAL(x), *end = xx + ncol * nrow, *f = REAL(F);
			for (; xx != end; xx += nrow, ++f)
				minimum<double>(xx, xx + nrow, *f);
			break;
		}
		default:
		{
			F = PROTECT(Rf_allocVector(INTSXP, ncol));
			int *xx = INTEGER(x), *end = xx + ncol * nrow, *f = INTEGER(F);
			for (; xx != end; xx += nrow, ++f)
				minimum<int>(xx, xx + nrow, *f);
			break;
		}
		}
		UNPROTECT(1);
		return F;
	}
}

// find the minimum index of its collumn
RcppExport SEXP Rfast_col_min_indices(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = col_min_indices(x);
	return __result;
	END_RCPP
}

// find the minimum value of its collumn
RcppExport SEXP Rfast_col_min(SEXP x, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	if (Rf_isMatrix(x))
		__result = col_min(x, parallel);
	else
		__result = Rfast::colMins(DataFrame(x), parallel);
	return __result;
	END_RCPP
}

// [[Rcpp::export]]
NumericVector row_min_indices(NumericMatrix x)
{
	const unsigned int p = x.nrow();
	mat X = mat(x.begin(), p, x.ncol(), false);
	NumericVector F(p);
	NumericVector::iterator FF = F.begin();
	for (unsigned int i = 0; i < p; ++i, ++FF)
		*FF = (X.row(i)).index_min() + 1;
	return F;
}

SEXP row_min(SEXP x)
{
	int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	SEXP F;
	F = PROTECT(Rf_allocVector(REALSXP, nrow));
	double *xx = REAL(x), *end = xx + ncol * nrow, *f = REAL(F), *x3, *ff;
	const double *endf = f + LENGTH(F);
	for (ff = f; ff != endf; ++ff, ++xx)
		*ff = *xx;
	for (; xx != end;)
		for (ff = f, x3 = xx, xx += nrow; x3 != xx; ++ff, ++x3)
		{
			*ff = std::min(*ff, *x3);
		}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_row_min_indices(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = row_min_indices(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_row_min(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = row_min(x);
	return __result;
	END_RCPP
}

////////////////////////////////////////////////////////////

SEXP col_min_max(SEXP x)
{
	const int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	SEXP F;
	switch (TYPEOF(x))
	{
	case REALSXP:
	{
		F = PROTECT(Rf_allocMatrix(REALSXP, 2, ncol));
		double *xx = REAL(x), *end = xx + LENGTH(x), *f = REAL(F), min, max;
		for (; xx != end; xx += nrow, f += 2)
		{
			min_max<double>(xx, xx + nrow, min, max);
			*f = min;
			f[1] = max;
		}
		break;
	}
	default:
	{
		F = PROTECT(Rf_allocMatrix(INTSXP, 2, ncol));
		int *xx = INTEGER(x), *end = xx + LENGTH(x), *f = INTEGER(F), min, max;
		for (; xx != end; xx += nrow, f += 2)
		{
			min_max<int>(xx, xx + nrow, min, max);
			*f = min;
			f[1] = max;
		}
	}
	}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_col_min_max(SEXP x, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	if (Rf_isMatrix(x))
		__result = col_min_max(x);
	else
		__result = Rfast::colMinsMaxs(DataFrame(x), parallel);
	return __result;
	END_RCPP
}

SEXP row_min_max(SEXP x)
{
	int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	SEXP F;
	F = PROTECT(Rf_allocMatrix(REALSXP, 2, nrow));
	double *xx = REAL(x), *end = xx + ncol * nrow, *f = REAL(F), *x3, *ff;
	const double *endf = f + (nrow << 1);
	for (ff = f; ff != endf; ff += 2, ++xx)
		*ff = ff[1] = *xx;
	for (; xx != end;)
		for (ff = f, x3 = xx, xx += nrow; x3 != xx; ff += 2, ++x3)
		{
			if (*ff > *x3)
				*ff = *x3;
			else if (ff[1] < *x3)
				ff[1] = *x3;
		}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_row_min_max(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = row_min_max(x);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////

using std::nth_element;

SEXP col_nth(NumericMatrix X, IntegerVector elems, const int num_of_nths, const bool descend, const bool na_rm, const bool index)
{
	const unsigned int n = elems.size();
	SEXP F = R_NilValue;
	if (num_of_nths == 1)
	{
		NumericVector y(X.nrow());
		if (index)
		{
			F = PROTECT(Rf_allocVector(INTSXP, n));
			int *ff = INTEGER(F);
			for (unsigned int i = 0; i != n; ++ff, ++i)
			{
				y = X.column(i);
				*ff = nth_helper_index<NumericVector>(y, elems[i], descend, na_rm);
			}
		}
		else
		{
			F = PROTECT(Rf_allocVector(REALSXP, n));
			double *ff = REAL(F);
			for (unsigned int i = 0; i != n; ++ff, ++i)
			{
				y = X.column(i);
				*ff = nth_helper<NumericVector>(y, elems[i], descend, na_rm);
			}
		}
	}
	else if (num_of_nths > 1)
	{
		colvec y(X.nrow());
		if (index)
		{
			F = PROTECT(Rf_allocMatrix(INTSXP, num_of_nths, n));
			NumericMatrix ff(F);
			mat x(X.begin(), X.nrow(), n, false), f(ff.begin(), num_of_nths, n, false);
			for (unsigned int i = 0; i < n; ++i)
			{
				y = x.col(i);
				f.col(i) = nth_helper_index_n_elems<colvec>(y, elems[i], descend, na_rm);
			}
		}
		else
		{
			F = PROTECT(Rf_allocMatrix(REALSXP, num_of_nths, n));
			NumericMatrix ff(F);
			mat x(X.begin(), X.nrow(), n, false), f(ff.begin(), num_of_nths, n, false);
			for (unsigned int i = 0; i < n; ++i)
			{
				y = x.col(i);
				f.col(i) = nth_helper_n_elems<colvec>(y, elems[i], descend, na_rm);
			}
		}
	}
	UNPROTECT(1);
	return F;
}

// nth_element
RcppExport SEXP Rfast_col_nth(SEXP xSEXP, SEXP ySEXP, SEXP num_of_nthsSEXP, SEXP descendSEXP, SEXP na_rmSEXP, SEXP indexSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<IntegerVector>::type y(ySEXP);
	traits::input_parameter<const unsigned int>::type num_of_nths(num_of_nthsSEXP);
	traits::input_parameter<const bool>::type descend(descendSEXP);
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	traits::input_parameter<const bool>::type index(indexSEXP);
	__result = col_nth(x, y, num_of_nths, descend, na_rm, index);
	return __result;
	END_RCPP
}

SEXP row_nth(NumericMatrix X, IntegerVector elems, const int num_of_nths, const bool descend, const bool na_rm, const bool index)
{
	const unsigned int n = elems.size();
	SEXP F = R_NilValue;

	if (num_of_nths == 1)
	{
		NumericVector y(X.ncol());
		if (index)
		{
			F = PROTECT(Rf_allocVector(INTSXP, n));
			int *ff = INTEGER(F);
			for (unsigned int i = 0; i != n; ++ff, ++i)
			{
				y = X.row(i);
				*ff = nth_helper_index<NumericVector>(y, elems[i], descend, na_rm);
			}
		}
		else
		{
			F = PROTECT(Rf_allocVector(REALSXP, n));
			double *ff = REAL(F);
			for (unsigned int i = 0; i != n; ++ff, ++i)
			{
				y = X.row(i);
				*ff = nth_helper<NumericVector>(y, elems[i], descend, na_rm);
			}
		}
	}
	else if (num_of_nths > 1)
	{
		rowvec y(X.ncol());
		if (index)
		{
			F = PROTECT(Rf_allocMatrix(INTSXP, n, num_of_nths));
			NumericMatrix ff(F);
			mat x(X.begin(), X.nrow(), n, false), f(ff.begin(), num_of_nths, n, false);
			for (unsigned int i = 0; i < n; ++i)
			{
				y = x.row(i);
				f.row(i) = nth_helper_index_n_elems<rowvec>(y, elems[i], descend, na_rm);
			}
		}
		else
		{
			F = PROTECT(Rf_allocMatrix(REALSXP, n, num_of_nths));
			NumericMatrix ff(F);
			mat x(X.begin(), X.nrow(), n, false), f(ff.begin(), num_of_nths, n, false);
			for (unsigned int i = 0; i < n; ++i)
			{
				y = x.row(i);
				f.row(i) = nth_helper_n_elems<rowvec>(y, elems[i], descend, na_rm);
			}
		}
	}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_row_nth(SEXP xSEXP, SEXP ySEXP, SEXP num_of_nthsSEXP, SEXP descendSEXP, SEXP na_rmSEXP, SEXP indexSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<IntegerVector>::type y(ySEXP);
	traits::input_parameter<const unsigned int>::type num_of_nths(num_of_nthsSEXP);
	traits::input_parameter<const bool>::type descend(descendSEXP);
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	traits::input_parameter<const bool>::type index(indexSEXP);
	__result = row_nth(x, y, num_of_nths, descend, na_rm, index);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////

IntegerMatrix col_order(NumericMatrix x, const bool stable, const bool descending)
{
	const int ncl = x.ncol();
	IntegerMatrix f(x.nrow(), ncl);
	for (int i = 0; i < ncl; ++i)
	{
		f.column(i) = Order(x.column(i), stable, descending,false);
	}
	return f;
}

RcppExport SEXP Rfast_col_order(SEXP xSEXP, SEXP stableSEXP, SEXP descendingSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const bool>::type stable(stableSEXP);
	traits::input_parameter<const bool>::type descending(descendingSEXP);
	__result = col_order(x, stable, descending);
	return __result;
	END_RCPP
}

IntegerMatrix row_order(NumericMatrix x, const bool stable, const bool descending)
{
	const int nrw = x.nrow();
	IntegerMatrix f(nrw, x.ncol());
	for (int i = 0; i < nrw; ++i)
		f.row(i) = Order(x.row(i), stable, descending,false);
	return f;
}

RcppExport SEXP Rfast_row_order(SEXP xSEXP, SEXP stableSEXP, SEXP descendingSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const bool>::type stable(stableSEXP);
	traits::input_parameter<const bool>::type descending(descendingSEXP);
	__result = row_order(x, stable, descending);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////

NumericVector col_prods(SEXP x, string method)
{
	const int n = Rf_ncols(x);
	NumericVector f(n);
	if (method == "direct")
	{
		mat X(REAL(x), Rf_nrows(x), n, false);
		rowvec ff(f.begin(), n, false);
		ff = prod(X, 0);
	}
	else if (method == "expsumlog")
	{
		int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
		double *xx = REAL(x), *end = xx + ncol * nrow, *ff = f.begin();
		for (; xx != end; xx += nrow, ++ff)
		{
			*ff = exp(accumulate(xx, xx + nrow, 0.0, [](double &s, double x)
								 { return x < 0 ? s + x : s + log(x); }));
		}
	}
	else
	{
		stop("Error: Unsupported method.");
	}
	return f;
}

RcppExport SEXP Rfast_col_prods(SEXP x, SEXP methodSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<string>::type method(methodSEXP);
	__result = col_prods(x, method);
	return __result;
	END_RCPP
}

NumericVector row_prods(NumericMatrix x)
{
	const int n = x.nrow();
	NumericVector f(n);
	mat X = mat(x.begin(), n, x.ncol(), false);
	colvec ff(f.begin(), n, false);
	ff = prod(X, 1);
	return f;
}

RcppExport SEXP Rfast_row_prods(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = row_prods(x);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////

using std::string;

DataFrame col_ranks(DataFrame x, string method, const bool descend, const bool stable, const bool parallel)
{
	List f(x.size());
	if (parallel)
	{
		if (method == "average")
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				colvec y;
				int i;
#pragma omp critical
				{
					NumericVector yy;
					yy = *s;
					y = colvec(yy.begin(), yy.size());
					i = s - x.begin();
				}
				y = rank_mean<colvec, colvec, ivec>(y, descend);
#pragma omp critical
				{
					f.insert(i, NumericVector(y.begin(), y.end()));
				}
			}
		}
		else if (method == "min")
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				colvec y;
				int i;
#pragma omp critical
				{
					NumericVector yy;
					yy = *s;
					y = colvec(yy.begin(), yy.size());
					i = s - x.begin();
				}
				y = rank_min<colvec, colvec, ivec>(y, descend);
#pragma omp critical
				{
					f.insert(i, NumericVector(y.begin(), y.end()));
				}
			}
		}
		else if (method == "max")
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				colvec y;
				int i;
#pragma omp critical
				{
					NumericVector yy;
					yy = *s;
					y = colvec(yy.begin(), yy.size());
					i = s - x.begin();
				}
				y = rank_max<colvec, colvec, ivec>(y, descend);
#pragma omp critical
				{
					f.insert(i, NumericVector(y.begin(), y.end()));
				}
			}
		}
		else if (method == "first")
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				colvec y;
				int i;
#pragma omp critical
				{
					NumericVector yy;
					yy = *s;
					y = colvec(yy.begin(), yy.size());
					i = s - x.begin();
				}
				y = rank_first<colvec, colvec, ivec>(y, descend, stable);
#pragma omp critical
				{
					f.insert(i, NumericVector(y.begin(), y.end()));
				}
			}
		}
		else
			stop("Error. Wrong method.");
	}
	else
	{
		int i = 0;
		NumericVector y(x.nrows());
		for (auto c : x)
		{
			y = c;
			f[i++] = Rank(y, method, descend, stable,false);
		}
	}
	f.names() = x.names();
	return f;
}

NumericMatrix col_ranks(NumericMatrix x, string method, const bool descend, const bool stable, const bool parallel)
{
	const int ncl = x.ncol(), nrw = x.nrow();
	NumericMatrix f(nrw, ncl);
	if (parallel)
	{
		mat xx(x.begin(), nrw, ncl, false);
		mat ff(f.begin(), nrw, ncl, false);
		if (method == "average")
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (int i = 0; i < ncl; ++i)
			{
				ff.col(i) = rank_mean<colvec, colvec, ivec>(xx.col(i), descend);
			}
		}
		else if (method == "min")
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (int i = 0; i < ncl; ++i)
			{
				ff.col(i) = rank_min<colvec, colvec, ivec>(xx.col(i), descend);
			}
		}
		else if (method == "max")
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (int i = 0; i < ncl; ++i)
			{
				ff.col(i) = rank_max<colvec, colvec, ivec>(xx.col(i), descend);
			}
		}
		else if (method == "first")
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (int i = 0; i < ncl; ++i)
			{
				ff.col(i) = rank_first<colvec, colvec, ivec>(xx.col(i), descend, stable);
			}
		}
		else
			stop("Error. Wrong method.");
	}
	else
	{
		for (int i = 0; i < ncl; ++i)
		{
			f.column(i) = Rank(x.column(i), method, descend, stable,false);
		}
	}
	return f;
}

RcppExport SEXP Rfast_col_ranks(SEXP xSEXP, SEXP methodSEXP, SEXP descendSEXP, SEXP stableSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type descend(descendSEXP);
	traits::input_parameter<const bool>::type stable(stableSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	if(Rf_isMatrix(xSEXP))
		__result = col_ranks(NumericMatrix(xSEXP), method, descend, stable, parallel);
	else
		__result = col_ranks(DataFrame(xSEXP), method, descend, stable, parallel);
	return __result;
	END_RCPP
}

NumericMatrix row_ranks(NumericMatrix x, string method, const bool descend, const bool stable)
{
	const int n = x.nrow();
	NumericMatrix f(n, x.ncol());
	for (int i = 0; i < n; ++i)
	{
		f.row(i) = Rank(x.row(i), method, descend, stable,false);
	}
	return f;
}

RcppExport SEXP Rfast_row_ranks(SEXP xSEXP, SEXP methodSEXP, SEXP descendSEXP, SEXP stableSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type descend(descendSEXP);
	traits::input_parameter<const bool>::type stable(stableSEXP);
	__result = row_ranks(x, method, descend, stable);
	return __result;
	END_RCPP
}

////////////////////////////////////////////////////////////////////////////

RcppExport SEXP Rfast_row_shuffle(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = Rfast::rowShuffle(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_shuffle(SEXP xSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	if(Rf_isMatrix(xSEXP))
		__result = Rfast::colShuffle(NumericMatrix(xSEXP));
	else
		__result = Rfast::colShuffle(DataFrame(xSEXP));
	return __result;
	END_RCPP
}

////////////////////////////////////////////////////////////////

template <class T, class Ret, class T1, class F1, class F2>
Ret col_sums(T1 x, SEXP indices, const bool na_rm = false)
{
	const int n = Rf_isNull(indices) ? 0 : LENGTH(indices);
	F1 X(x.begin(), x.nrow(), x.ncol(), false);
	Ret f(n == 0 ? X.n_cols : n);
	if (n == 0)
	{
		F2 ff(f.begin(), X.n_cols, false);
		if (na_rm)
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (unsigned int i = 0; i < X.n_cols; ++i)
			{
				ff[i] = sum_with_condition<T, notNA<T>, typename F1::col_iterator>(X.begin_col(i), X.end_col(i));
			}
		}
		else
		{
			ff = sum(X, 0);
		}
	}
	else
	{
		IntegerVector ind(indices);
		if (na_rm)
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (unsigned int i = 0; i < X.n_cols; ++i)
			{
				const int j = ind[i] - 1;
				f[i] = sum_with_condition<T, notNA<T>, typename F1::col_iterator>(X.begin_col(j), X.end_col(j));
			}
		}
		else
		{
			for (int i = 0; i < n; ++i)
			{
				const int j = ind[i] - 1;
				f[i] = accumulate(X.begin_col(j), X.end_col(j), 0);
			}
		}
	}
	return f;
}

template <class T, class Ret, class T1, class F1, class F2>
Ret row_sums(T1 x, SEXP indices, const bool na_rm = false)
{
	const int n = Rf_isNull(indices) ? 0 : LENGTH(indices);
	F1 X(x.begin(), x.nrow(), x.ncol(), false);
	Ret f(n == 0 ? X.n_rows : n);
	if (n == 0)
	{
		F2 ff(f.begin(), X.n_rows, false);
		if (na_rm)
		{
#pragma omp parallel for
			for (unsigned int i = 0; i < X.n_rows; ++i)
			{
				ff[i] = sum_with_condition<T, notNA<T>, typename F1::row_iterator>(X.begin_row(i), X.end_row(i));
			}
		}
		else
		{
			ff = sum(X, 1);
		}
	}
	else
	{
		IntegerVector ind(indices);
		if (na_rm)
		{
#pragma omp parallel for
			for (unsigned int i = 0; i < X.n_rows; ++i)
			{
				const int j = ind[i] - 1;
				f[i] = sum_with_condition<T, notNA<T>, typename F1::row_iterator>(X.begin_row(j), X.end_row(j));
			}
		}
		else
		{
			for (int i = 0; i < n; ++i)
			{
				const int j = ind[i] - 1;
				f[i] = accumulate(X.begin_col(j), X.end_col(j), 0);
			}
		}
	}
	return f;
}

RcppExport SEXP Rfast_row_sums(SEXP x, SEXP indices, SEXP na_rmSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	__result = Rf_isInteger(x) ? row_sums<int, IntegerVector, IntegerMatrix, imat, icolvec>(x, indices, na_rm)
							   : row_sums<double, NumericVector, NumericMatrix, mat, colvec>(x, indices, na_rm);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_sums(SEXP x, SEXP indices, SEXP na_rmSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	__result = Rf_isInteger(x) ? col_sums<int, IntegerVector, IntegerMatrix, imat, irowvec>(x, indices, na_rm)
							   : col_sums<double, NumericVector, NumericMatrix, mat, rowvec>(x, indices, na_rm);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////

IntegerMatrix col_tabulate(IntegerMatrix x, int nroww)
{
	const int ncl = x.ncol();
	IntegerMatrix f(nroww, ncl);
	for (int i = 0; i < ncl; ++i)
		f.column(i) = Tabulate<IntegerVector>(x.column(i), nroww);
	return f;
}

IntegerMatrix row_tabulate(IntegerMatrix x, int ncoll)
{
	const int nrw = x.nrow();
	IntegerMatrix f(nrw, ncoll);
	for (int i = 0; i < nrw; ++i)
		f.row(i) = Tabulate<IntegerVector>(x.row(i), ncoll);
	return f;
}

RcppExport SEXP Rfast_row_tabulate(SEXP xSEXP, SEXP ncollSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<IntegerMatrix>::type x(xSEXP);
	traits::input_parameter<int>::type ncoll(ncollSEXP);
	__result = row_tabulate(x, ncoll);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_tabulate(SEXP xSEXP, SEXP nrowwSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<IntegerMatrix>::type x(xSEXP);
	traits::input_parameter<int>::type nroww(nrowwSEXP);
	__result = col_tabulate(x, nroww);
	return __result;
	END_RCPP
}

////////////////////////////////////////////////////////////

SEXP col_true(SEXP x)
{
	const int p = Rf_nrows(x);
	SEXP f = Rf_allocVector(INTSXP, p);
	int *ff = INTEGER(f), *xx = LOGICAL(x), *endx = xx + LENGTH(x);
	for (; xx != endx; xx += p, ++ff)
		*ff = True(xx, xx + p);
	return f;
}

SEXP row_true(SEXP x)
{
	int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	SEXP F = PROTECT(Rf_allocVector(INTSXP, nrow));
	int *xx = INTEGER(x), *end = xx + ncol * nrow, *f = INTEGER(F), *startx, *startf;
	const int *endf = f + LENGTH(F);
	for (startf = f; startf != endf; ++startf)
		*startf = 0;
	while (xx != end)
	{
		for (startf = f, startx = xx, xx += nrow; startx != xx; ++startf, ++startx)
		{
			*startf += *startx;
		}
	}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_row_true(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = row_true(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_true(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_true(x);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////////

SEXP col_true_false(SEXP x)
{
	const int n = Rf_nrows(x);
	SEXP f = Rf_allocMatrix(INTSXP, 2, Rf_ncols(x));
	int *ff = INTEGER(f), *xx = LOGICAL(x), *endx = xx + LENGTH(x), t;
	for (; xx != endx; xx += n, ff += 2)
	{
		t = True(xx, xx + n);
		*ff = n - t;
		ff[1] = t;
	}
	return f;
}

RcppExport SEXP Rfast_col_true_false(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_true_false(x);
	return __result;
	END_RCPP
}

SEXP row_true_false(SEXP x)
{
	int ncol = Rf_ncols(x), nrow = Rf_nrows(x);
	SEXP F = PROTECT(Rf_allocMatrix(INTSXP, 2, nrow));
	int *xx = INTEGER(x), *end = xx + ncol * nrow, *f = INTEGER(F), *startx, *startf;
	const int *endf = f + LENGTH(F);
	for (startf = f; startf != endf; startf += 2)
	{
		*startf = ncol;
		startf[1] = 0;
	}
	while (xx != end)
	{
		for (startf = f, startx = xx, xx += nrow; startx != xx; startf += 2, ++startx)
		{
			*startf -= *startx;
			startf[1] += *startx;
		}
	}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_row_true_false(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = row_true_false(x);
	return __result;
	END_RCPP
}

////////////////////////////////////////////////

SEXP col_pmin(SEXP x, SEXP y)
{
	const int nrows = Rf_nrows(x), ncols = Rf_ncols(x);
	SEXP f = Rf_allocMatrix(REALSXP, nrows, ncols);
	double *startx = REAL(x), *end = startx + ncols * nrows, *starty = REAL(y), *startf = REAL(f), *endx;
	for (; startx != end;)
	{
		endx = startx + nrows;
		for (; startx != endx; ++startx, ++starty, ++startf)
			*startf = std::min(*startx, *starty);
	}
	return f;
}

RcppExport SEXP Rfast_col_pmin(SEXP x, SEXP y)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_pmin(x, y);
	return __result;
	END_RCPP
}

//////////////////////////////////////////////////////////////

SEXP col_pmax(SEXP x, SEXP y)
{
	const int nrows = Rf_nrows(x), ncols = Rf_ncols(x);
	SEXP f = Rf_allocMatrix(REALSXP, nrows, ncols);
	double *startx = REAL(x), *end = startx + ncols * nrows, *starty = REAL(y), *startf = REAL(f), *endx;
	for (; startx != end;)
	{
		endx = startx + nrows;
		for (; startx != endx; ++startx, ++starty, ++startf)
			*startf = std::max(*startx, *starty);
	}
	return f;
}

RcppExport SEXP Rfast_col_pmax(SEXP x, SEXP y)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_pmax(x, y);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////

SEXP col_cum_prods(SEXP x)
{
	const int p = Rf_nrows(x);
	SEXP f = Rf_duplicate(x);
	double *ff = REAL(f), *endf = ff + LENGTH(f);
	int i = 1;
	for (++ff; ff != endf; ++ff, ++i)
	{
		i != p ? *ff *= ff[-1] : i = 0;
	}
	return f;
}

RcppExport SEXP Rfast_col_cum_prods(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_cum_prods(x);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////

SEXP col_cum_maxs(SEXP x)
{
	int p = Rf_nrows(x);
	SEXP f = Rf_duplicate(x);
	double *ff = REAL(f), *endf = ff + LENGTH(f);
	int i = 1;
	for (++ff; ff != endf; ++ff, ++i)
	{
		i != p ? *ff = std::max(*ff, ff[-1]) : i = 0;
	}
	return f;
}

RcppExport SEXP Rfast_col_cum_maxs(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_cum_maxs(x);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////

SEXP col_cum_sums(SEXP x)
{
	int p = Rf_nrows(x);
	SEXP f = Rf_duplicate(x);
	double *ff = REAL(f), *endf = ff + LENGTH(f);
	int i = 1;
	for (++ff; ff != endf; ++ff, ++i)
	{
		i != p ? *ff += ff[-1] : i = 0;
	}
	return f;
}

RcppExport SEXP Rfast_col_cum_sums(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_cum_sums(x);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////

SEXP col_cum_mins(SEXP x)
{
	int p = Rf_nrows(x);
	SEXP f = Rf_duplicate(x);
	double *ff = REAL(f), *endf = ff + LENGTH(f);
	int i = 1;
	for (++ff; ff != endf; ++ff, ++i)
	{
		i != p ? *ff = std::min(*ff, ff[-1]) : i = 0;
	}
	return f;
}

RcppExport SEXP Rfast_col_cum_mins(SEXP x)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = col_cum_mins(x);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////

bool col_row_value(NumericMatrix X, double v)
{
	int i, n = X.nrow(), p = X.ncol();
	mat x(X.begin(), n, p, false);
	for (i = 0; i < p; ++i)
	{
		if (all(x.col(i) == v))
		{
			return true;
		}
	}
	for (i = 0; i < n; ++i)
	{
		if (all(x.row(i) == v))
		{
			return true;
		}
	}
	return false;
}

RcppExport SEXP Rfast_col_row_value(SEXP xSEXP, SEXP vSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const double>::type v(vSEXP);
	__result = col_row_value(x, v);
	return __result;
	END_RCPP
}

//////////////////////////////////////////////////////////////////////////////////

NumericMatrix columns(NumericMatrix x, IntegerVector ind)
{
	const int nrw = x.nrow(), ncl = ind.size();
	NumericMatrix f(nrw, ncl);
	mat ff(f.begin(), nrw, ncl, false), xx(x.begin(), nrw, x.ncol(), false);
	for (int i = 0; i < ncl; ++i)
		ff.col(i) = xx.col(ind[i] - 1);
	return f;
}

SEXP rows(SEXP X, SEXP Ind)
{
	const int nrw = Rf_nrows(X), ncl = Rf_ncols(X);
	SEXP F = PROTECT(Rf_allocMatrix(REALSXP, LENGTH(Ind), ncl));
	double *start = REAL(X), *ff = REAL(F), *xx = start;
	int *start_ind = INTEGER(Ind), *ind, *end_ind = start_ind + LENGTH(Ind);
	for (int i = 0; i < ncl; ++i)
	{
		for (ind = start_ind; ind != end_ind; ++ind)
		{
			xx = start + *ind - 1;
			*ff++ = *xx;
		}
		start += nrw;
	}
	UNPROTECT(1);
	return F;
}

RcppExport SEXP Rfast_columns(SEXP xSEXP, SEXP indSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<IntegerVector>::type ind(indSEXP);
	__result = columns(x, ind);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_rows(SEXP x, SEXP ind)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	__result = rows(x, ind);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_col_vars(SEXP xSEXP, SEXP stdSEXP, SEXP na_rmSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const bool>::type std(stdSEXP);
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = Rf_isMatrix(xSEXP) ? Rfast::colVars(NumericMatrix(xSEXP), std, na_rm, parallel) : Rfast::colVars(DataFrame(xSEXP), std, na_rm, parallel);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_row_vars(SEXP xSEXP, SEXP stdSEXP, SEXP na_rmSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const bool>::type std(stdSEXP);
	traits::input_parameter<const bool>::type na_rm(na_rmSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = Rfast::rowVars(x, std, na_rm, parallel);
	return __result;
	END_RCPP
}
