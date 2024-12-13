
// Author: Manos Papadakis
//[[Rcpp::plugins(cpp11)]]

#define ARMA_64BIT_WORD

#include <RcppArmadillo.h>

#include "Dist.h"
#include "mn.h"

using namespace Rcpp;
using namespace arma;

List dcor(NumericMatrix &X, NumericMatrix &Y, const long double &n2, const long double &n3, const long double &n4) {
	const size_t ncl = X.ncol(), nrw = X.nrow();
	mat x(X.begin(), nrw, ncl, false), y(Y.begin(), nrw, ncl, false);
	colvec sum_row_sa(ncl), sum_row_sb(ncl);
	long double sum_sa = 0.0, sum_sb = 0.0, sum_sab = 0.0, sum_sa2 = 0.0, sum_sb2 = 0.0;
	for (size_t i = 0; i < ncl - 1; ++i) {
		colvec xv(x.begin_col(i), nrw, false), yv(y.begin_col(i), nrw, false);
		long double a = 0.0, b = 0.0;
		for (size_t j = i + 1; j < ncl; ++j) {
			colvec xj(x.begin_col(j), nrw, false), yj(y.begin_col(j), nrw, false);
			long double vx = Dist::euclidean<false>(xv, xj);
			long double vy = Dist::euclidean<false>(yv, yj);
			a += vx;
			b += vy;
			sum_sa += vx;
			sum_sb += vy;
			sum_sab += vx * vy;
			sum_sa2 += vx * vx;
			sum_sb2 += vy * vy;
			sum_row_sa[j] += vx;
			sum_row_sb[j] += vy;
		}
		sum_row_sa[i] += a;
		sum_row_sb[i] += b;
	}
	sum_sa *= 2;
	sum_sb *= 2;
	sum_sa2 *= 2;
	sum_sb2 *= 2;
	sum_sab *= 2;

	long double dvarX = sqrt(sum_sa2 / n2 - 2 * dot(sum_row_sa, sum_row_sa) / n3 + sum_sa * sum_sa / n4);
	long double dvarY = sqrt(sum_sb2 / n2 - 2 * dot(sum_row_sb, sum_row_sb) / n3 + sum_sb * sum_sb / n4);
	long double dcov = sqrt(sum_sab / n2 - 2 * dot(sum_row_sa, sum_row_sb) / n3 + sum_sa * sum_sb / n4);
	long double dcor = dcov / sqrt(dvarX * dvarY);

	return List::create(_["dcov"] = dcov, _["dvarX"] = dvarX, _["dvarY"] = dvarY, _["dcor"] = dcor);
}


List dcor(NumericMatrix X, NumericMatrix Y) {
	const double nd = X.ncol();
	const long double n2 = nd * nd, n3 = n2 * nd, n4 = n3 * nd;
	return dcor(X, Y, n2, n3, n4);
}

List bcdcor(NumericMatrix X, NumericMatrix Y) {
	const double nd = X.ncol();
	const long double n2 = nd * (nd - 3), n3 = n2 * (nd - 2), n4 = n3 * (nd - 1);
	return dcor(X, Y, n2, n3, n4);
}

RcppExport SEXP Rfast_dcor(SEXP xSEXP, SEXP ySEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericMatrix>::type y(ySEXP);
	__result = dcor(x, y);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_bcdcor(SEXP xSEXP, SEXP ySEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericMatrix>::type y(ySEXP);
	__result = bcdcor(x, y);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////////

double dcov(NumericMatrix &X, NumericMatrix &Y, const long double &n2, const long double &n3, const long double &n4) {
	const size_t ncl = X.ncol(), nrw = X.nrow();
	mat x(X.begin(), nrw, ncl, false), y(Y.begin(), nrw, ncl, false);
	colvec sum_row_sa(ncl), sum_row_sb(ncl);
	long double sum_sa = 0.0, sum_sb = 0.0, sum_sab = 0.0;
	for (size_t i = 0; i < ncl - 1; ++i) {
		colvec xv(x.begin_col(i), nrw, false), yv(y.begin_col(i), nrw, false);
		long double a = 0.0, b = 0.0;
		for (size_t j = i + 1; j < ncl; ++j) {
			colvec xj(x.begin_col(j), nrw, false), yj(y.begin_col(j), nrw, false);
			long double vx = Dist::euclidean<false>(xv, xj);
			long double vy = Dist::euclidean<false>(yv, yj);
			sum_sa += vx;
			sum_sb += vy;
			a += vx;
			b += vy;
			sum_sab += vx * vy;
			sum_row_sa[j] += vx;
			sum_row_sb[j] += vy;
		}
		sum_row_sa[i] += a;
		sum_row_sb[i] += b;
	}
	sum_sa *= 2;
	sum_sb *= 2;
	sum_sab *= 2;
	return sqrt(sum_sab / n2 - 2 * dot(sum_row_sa, sum_row_sb) / n3 + sum_sa * sum_sb / n4);
}

double dcov(NumericMatrix X, NumericMatrix Y) {
	const double nd = X.ncol();
	const long double n2 = nd * nd, n3 = n2 * nd, n4 = n3 * nd;
	return dcov(X, Y, n2, n3, n4);
}

double bcdcov(NumericMatrix X, NumericMatrix Y) {
	const double nd = X.ncol();
	const long double n2 = nd * (nd - 3), n3 = n2 * (nd - 2), n4 = n3 * (nd - 1);
	return dcov(X, Y, n2, n3, n4);
}

RcppExport SEXP Rfast_dcov(SEXP xSEXP, SEXP ySEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericMatrix>::type y(ySEXP);
	__result = dcov(x, y);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_bcdcov(SEXP xSEXP, SEXP ySEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericMatrix>::type y(ySEXP);
	__result = bcdcov(x, y);
	return __result;
	END_RCPP
}

////////////////////////////////////////////////////////////////////

double dvar(NumericMatrix &X, const long double &n2, const long double &n3, const long double &n4) {
	const int n = X.ncol();
	const size_t ncl = X.ncol(), nrw = X.nrow();
	mat x(X.begin(), nrw, ncl, false);
	colvec sum_row_sa(n);
	long double sum_sa = 0.0, sum_sa2 = 0.0;
	for (size_t i = 0; i < ncl - 1; ++i) {
		colvec xv(x.begin_col(i), nrw, false);
		long double a = 0.0;
		for (size_t j = i + 1; j < ncl; ++j) {
			colvec y(x.begin_col(j), nrw, false);
			long double v = Dist::euclidean<false>(xv, y);
			sum_sa += v;
			a += v;
			sum_row_sa[j] += v;
			sum_sa2 += v * v;
		}
		sum_row_sa[i] += a;
	}
	sum_sa2 *= 2;
	sum_sa *= 2;
	return sqrt(sum_sa2 / n2 - 2 * sum(square(sum_row_sa)) / n3 + sum_sa * sum_sa / n4);
}

double dvar(NumericMatrix X) {
	const double nd = X.ncol();
	const long double n2 = nd * nd, n3 = n2 * nd, n4 = n3 * nd;
	return dvar(X, n2, n3, n4);
}


double bcdvar(NumericMatrix X) {
	const double nd = X.ncol();
	const long double n2 = nd * (nd - 3), n3 = n2 * (nd - 2), n4 = n3 * (nd - 1);
	return dvar(X, n2, n3, n4);
}

RcppExport SEXP Rfast_dvar(SEXP xSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = dvar(x);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_bcdvar(SEXP xSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	__result = bcdvar(x);
	return __result;
	END_RCPP
}

//////////////////////////////////////////////////////////////////

double edist(NumericMatrix x, NumericMatrix y) {
	const int n1 = x.ncol(), n2 = y.ncol();
	double mij = total_dista(x, y, "euclidean", false), mii = total_dist(x, "euclidean", false),
		   mjj = total_dist(y, "euclidean", false);
	return (2 * mij / (n1 * n2) - 2 * mii / (n1 * n1) - 2 * mjj / (n2 * n2)) * n1 * n2 / (n1 + n2);
}

RcppExport SEXP Rfast_edist(SEXP xSEXP, SEXP ySEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericMatrix>::type y(ySEXP);
	__result = edist(x, y);
	return __result;
	END_RCPP
}

#undef ARMA_64BIT_WORD

/////////////////////////////////////////////////////////////////////////////////////
