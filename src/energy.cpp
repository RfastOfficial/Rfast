
// Author: Manos Papadakis
//[[Rcpp::plugins(cpp11)]]

#define ARMA_64BIT_WORD

#include <RcppArmadillo.h>

#include "Dist.h"
#include "mn.h"

using namespace Rcpp;
using namespace arma;

List dcor(NumericMatrix X, NumericMatrix Y) {
	const size_t ncl = X.ncol(), nrw = X.nrow();
	const double nd = ncl;
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

	long double n2 = nd * nd, n3 = n2 * nd, n4 = n3 * nd;

	long double dvarX = sqrt(sum_sa2 / n2 - 2 * dot(sum_row_sa, sum_row_sa) / n3 + sum_sa * sum_sa / n4);
	long double dvarY = sqrt(sum_sb2 / n2 - 2 * dot(sum_row_sb, sum_row_sb) / n3 + sum_sb * sum_sb / n4);
	long double dcov = sqrt(sum_sab / n2 - 2 * dot(sum_row_sa, sum_row_sb) / n3 + sum_sa * sum_sb / n4);
	long double dcor = dcov / sqrt(dvarX * dvarY);

	return List::create(_["dcov"] = dcov, _["dvarX"] = dvarX, _["dvarY"] = dvarY, _["dcor"] = dcor);
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

///////////////////////////////////////////////////////////////////////////

double dcov(NumericMatrix X, NumericMatrix Y) {
	const size_t ncl = X.ncol(), nrw = X.nrow();
	const double nd = ncl;
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
	return sqrt(sum_sab / (nd * nd) - 2 * dot(sum_row_sa, sum_row_sb) / (nd * nd * nd) +
				sum_sa * sum_sb / (nd * nd * nd * nd));
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

////////////////////////////////////////////////////////////////////

double dvar(NumericMatrix X) {
	const int n = X.ncol();
	const double nd = n;
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
	return sqrt(sum_sa2 / (nd * nd) - 2 * sum(square(sum_row_sa)) / (nd * nd * nd) +
				sum_sa * sum_sa / (nd * nd * nd * nd));
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

//////////////////////////////////////////////////////////////////

double edist(NumericMatrix x, NumericMatrix y) {
	const int n1 = x.ncol(), n2 = y.ncol();
	double mij = total_dista(x, y, "euclidean", false), mii = DistTotal::euclidean(x, false),
		   mjj = DistTotal::euclidean(y, false);
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

///////////////////////////////////////////////////////////////////////////////////////

static vec lower_tri(mat &x) {
	int n = x.n_cols, i, j;
	vec f(n * (n - 1) * 0.5);
	vec::iterator ff = f.begin();
	for (i = 0; i < n; ++i)
		for (j = i + 1; j < n; ++j, ++ff) *ff = x(j, i);
	return f;
}

double bcdcor(NumericMatrix x, NumericMatrix y) {
	NumericMatrix a = dist(x, "euclidean", false);
	NumericMatrix b = dist(y, "euclidean", false);
	const double n = a.ncol();
	mat aa(a.begin(), n, n, false);
	mat bb(b.begin(), n, n, false);
	rowvec ma = mean(aa, 0), dgA(n);
	rowvec mb = mean(bb, 0), dgB(n);
	const double mean_ma = mean(ma), mean_mb = mean(mb);
	double n_1 = n / (n - 1), n_2 = n / (n - 2);
	dgA = ma - mean_ma;
	dgB = mb - mean_mb;
	mat A = aa.each_row() - ma;
	A = A.each_col() - ma.t();
	A = A + mean_ma - aa / n;
	vec Al = lower_tri(A);
	mat B = bb.each_row() - mb;
	B = B.each_col() - mb.t();
	B = B + mean_mb - bb / n;
	vec Bl = lower_tri(B);
	n_1 *= n_1;
	n_2 *= n_1;
	const double sdgab = sum(dgA % dgB), sdga = sum(square(dgA)), sdgb = sum(square(dgB));
	const double XY = n_1 * (sum(Al % Bl) * 2 + sdgab) - n_2 * sdgab;
	const double XX = n_1 * (sum(square(Al)) * 2 + sdga) - n_2 * sdga;
	const double YY = n_1 * (sum(square(Bl)) * 2 + sdgb) - n_2 * sdgb;
	return XY / sqrt(XX * YY);
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

#undef ARMA_64BIT_WORD

/////////////////////////////////////////////////////////////////////////////////////
