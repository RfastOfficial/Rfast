// Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "mn.h"

using namespace Rcpp;

//[[Rcpp::export]]
colvec rmdp(NumericMatrix Y, const int h, umat rnd, const int itertime, const bool parallel = false)
{
	const int n = Y.nrow();
	mat y(Y.begin(), n, Y.ncol(), false);
	double bestdet = 0;
	colvec final_vec(n);
	span index(0, h - 1);

	if (parallel)
	{
#ifdef _OPENMP
	#pragma omp parallel for
#endif
		for (int A = 0; A < itertime; ++A)
		{
			mat ny = y.rows(rnd(span::all, A));
			rowvec mu_t = mean(ny, 0);
			rowvec var_t = colvar_rmdp(ny);
			mat tmp = y.each_row() - mu_t;
			tmp = square(tmp);
			mat sama = tmp.each_row() / var_t;
			colvec disa = sum(sama, 1);
			colvec jvec(n, fill::zeros), ivec(n);
			uvec dist_perm(n), indextony(h), t(h);
			for (int l = 0, crit = 10; crit && l <= 15;)
			{
				l++;
				ivec.fill(0);
				dist_perm = Order_rmdp(disa);
				t = dist_perm(index);
				indextony(index) = t;
				ivec(t).fill(1);
				crit = accu(abs(ivec - jvec));
				jvec = ivec;
				ny = y.rows(indextony);
				mu_t = mean(ny, 0);
				var_t = var(ny, 0, 0);
				tmp = y.each_row() - mu_t;
				tmp = square(tmp);
				sama = tmp.each_row() / var_t;
				disa = sum(sama, 1);
			}

			double tempdet = prod(var_t);
			
#ifdef _OPENMP
#pragma omp critical
			{
#endif
				if (!bestdet || tempdet < bestdet)
				{
					bestdet = tempdet;
					final_vec = jvec;
				}
#ifdef _OPENMP
			}
#endif
		}
	}
	else
	{
		mat ny, tmp, sama;
		colvec jvec(n, fill::zeros), ivec(n), disa(n);
		uvec dist_perm(n), indextony(h), t(h);
		rowvec mu_t, var_t;
		double tempdet = 0.0;
		for (int A = 0; A < itertime; ++A)
		{
			ny = y.rows(rnd(span::all, A));
			mu_t = mean(ny, 0);
			var_t = colvar_rmdp(ny);
			tmp = y.each_row() - mu_t;
			tmp = square(tmp);
			sama = tmp.each_row() / var_t;
			disa = sum(sama, 1);
			for (int l = 0, crit = 10; crit && l <= 15;)
			{
				l++;
				ivec.fill(0);
				dist_perm = Order_rmdp(disa);
				t = dist_perm(index);
				indextony(index) = t;
				ivec(t).fill(1);
				crit = accu(abs(ivec - jvec));
				jvec = ivec;
				ny = y.rows(indextony);
				mu_t = mean(ny, 0);
				var_t = var(ny, 0, 0);
				tmp = y.each_row() - mu_t;
				tmp = square(tmp);
				sama = tmp.each_row() / var_t;
				disa = sum(sama, 1);
			}
			tempdet = prod(var_t);
			if (!bestdet || tempdet < bestdet)
			{
				bestdet = tempdet;
				final_vec = jvec;
			}
		}
	}
	return final_vec;
}

RcppExport SEXP Rfast_rmdp(SEXP ySEXP, SEXP hSEXP, SEXP rndSEXP, SEXP itertimeSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type y(ySEXP);
	traits::input_parameter<const int>::type h(hSEXP);
	traits::input_parameter<umat>::type rnd(rndSEXP);
	traits::input_parameter<const int>::type itertime(itertimeSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = rmdp(y, h, rnd, itertime, parallel);
	return __result;
	END_RCPP
}