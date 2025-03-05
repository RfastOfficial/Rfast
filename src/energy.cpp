
// Author: Manos Papadakis
//[[Rcpp::plugins(cpp11)]]

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

#include "Rfast/Dist.h"
#include "mn.h"

using namespace Rcpp;
using namespace arma;
using namespace Rfast;


void initSizes(const size_t n, long double &n2, long double &n3, long double &n4, const bool bc = false){
	const double nd = n;

	if(bc){
		n2 = nd * (nd - 3), n3 = n2 * (nd - 2), n4 = n3 * (nd - 1);
	}else{
		n2 = nd * nd, n3 = n2 * nd, n4 = n3 * nd;
	}
}


///////////////////////////////////////////////////////////////////////////

double dcov(const double &sum_sa, const double &sum_sb, const double &sum_sab, colvec &sum_row_sa, colvec &sum_row_sb, 
			const long double &n2, const long double &n3, const long double &n4, const bool bc = false) {
	auto res = sum_sab / n2 - 2 * dot(sum_row_sa, sum_row_sb) / n3 + sum_sa * sum_sb / n4;

	if(!bc){
		res = sqrt(res);
	}

	return res;
}

double dcov(NumericMatrix X, NumericMatrix Y, const bool bc = false) {
	const size_t ncl = X.ncol(), nrw = X.nrow();
	mat x(X.begin(), nrw, ncl, false), y(Y.begin(), nrw, ncl, false);
	colvec sum_row_sa(ncl), sum_row_sb(ncl);
	long double sum_sa = 0.0, sum_sb = 0.0, sum_sab = 0.0;
	for (size_t i = 0; i < ncl - 1; ++i) {
		colvec xv(x.begin_col(i), nrw, false), yv(y.begin_col(i), nrw, false);
		long double a = 0.0, b = 0.0;
		for (size_t j = i + 1; j < ncl; ++j) {
			colvec xj(x.begin_col(j), nrw, false), yj(y.begin_col(j), nrw, false);
			long double vx = Rfast::Dist::euclidean<true>(xv, xj);
			long double vy = Rfast::Dist::euclidean<true>(yv, yj);
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
	
	long double n2, n3, n4;

	initSizes(ncl, n2, n3, n4, bc);

	return dcov(sum_sa, sum_sb, sum_sab, sum_row_sa, sum_row_sb ,n2, n3, n4, bc);
}

RcppExport SEXP Rfast_dcov(SEXP xSEXP, SEXP ySEXP, SEXP bcSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericMatrix>::type y(ySEXP);
	traits::input_parameter<const bool>::type bc(bcSEXP);
	__result = dcov(x, y, bc);
	return __result;
	END_RCPP
}

////////////////////////////////////////////////////////////////////


double dvar(const double &sum_sa, const double &sum_sa2, colvec &sum_row_sa,
			const long double &n2, const long double &n3, const long double &n4, const bool bc = false) {
	auto res = sum_sa2 / n2 - 2 * sum(square(sum_row_sa)) / n3 + sum_sa * sum_sa / n4;
	
	if(!bc){
		res = sqrt(res);
	}

	return res;
}

double dvar(NumericMatrix X, const bool bc = false) {
	const size_t ncl = X.ncol(), nrw = X.nrow();
	mat x(X.begin(), nrw, ncl, false);
	colvec sum_row_sa(ncl);
	long double sum_sa = 0.0, sum_sa2 = 0.0;
	for (size_t i = 0; i < ncl - 1; ++i) {
		colvec xv(x.begin_col(i), nrw, false);
		long double a = 0.0;
		for (size_t j = i + 1; j < ncl; ++j) {
			colvec y(x.begin_col(j), nrw, false);
			long double v = Rfast::Dist::euclidean<true>(xv, y);
			sum_sa += v;
			a += v;
			sum_row_sa[j] += v;
			sum_sa2 += v * v;
		}
		sum_row_sa[i] += a;
	}
	sum_sa2 *= 2;
	sum_sa *= 2;
	
	long double n2, n3, n4;

	initSizes(ncl, n2, n3, n4, bc);

	return dvar(sum_sa, sum_sa2, sum_row_sa ,n2, n3, n4, bc);
}

RcppExport SEXP Rfast_dvar(SEXP xSEXP, SEXP bcSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const bool>::type bc(bcSEXP);
	__result = dvar(x, bc);
	return __result;
	END_RCPP
}

List dcor(const double &sum_sa, const double &sum_sa2, const double &sum_sb, const double &sum_sb2, const double &sum_sab, colvec &sum_row_sa, colvec &sum_row_sb, 
			const long double &n2, const long double &n3, const long double &n4, const bool bc = false){
	long double dvarX = dvar(sum_sa, sum_sa2, sum_row_sa ,n2, n3, n4, bc);
	long double dvarY = dvar(sum_sb, sum_sb2, sum_row_sb ,n2, n3, n4, bc);
	long double dcovv = dcov(sum_sa, sum_sb, sum_sab, sum_row_sa, sum_row_sb ,n2, n3, n4, bc);
	long double dcor = dcovv / sqrt(dvarX * dvarY);

	return List::create(_["dcov"] = dcovv, _["dvarX"] = dvarX, _["dvarY"] = dvarY, _["dcor"] = dcor);
}

List dcor(NumericMatrix X, NumericMatrix Y, const bool bc = false) {
	const size_t ncl = X.ncol(), nrw = X.nrow();
	mat x(X.begin(), nrw, ncl, false), y(Y.begin(), nrw, ncl, false);
	colvec sum_row_sa(ncl), sum_row_sb(ncl);
	long double sum_sa = 0.0, sum_sb = 0.0, sum_sab = 0.0, sum_sa2 = 0.0, sum_sb2 = 0.0;
	for (size_t i = 0; i < ncl - 1; ++i) {
		colvec xv(x.begin_col(i), nrw, false), yv(y.begin_col(i), nrw, false);
		long double a = 0.0, b = 0.0;
		for (size_t j = i + 1; j < ncl; ++j) {
			colvec xj(x.begin_col(j), nrw, false), yj(y.begin_col(j), nrw, false);
			long double vx = Rfast::Dist::euclidean<true>(xv, xj);
			long double vy = Rfast::Dist::euclidean<true>(yv, yj);
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

	long double n2, n3, n4;

	initSizes(ncl, n2, n3, n4, bc);

	return dcor(sum_sa, sum_sa2, sum_sb, sum_sb2, sum_sab, sum_row_sa, sum_row_sb ,n2, n3, n4, bc);
}


RcppExport SEXP Rfast_dcor(SEXP xSEXP, SEXP ySEXP, SEXP bcSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericMatrix>::type y(ySEXP);
	traits::input_parameter<const bool>::type bc(bcSEXP);
	__result = dcor(x, y, bc);
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

/////////////////////////////////////////////////////////////////////////////////////