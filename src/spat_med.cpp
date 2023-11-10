// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "Rfast.h"
#include "mn.h"

using namespace arma;
using namespace Rcpp;

rowvec spat_med(NumericMatrix Y, const double tol = 1e-09)
{
	mat y = mat(Y.begin(), Y.nrow(), Y.ncol(), false);
	mat z(y.n_rows, y.n_cols);
	colvec ww(z.n_rows);
	rowvec u1(y.n_rows), u2;
	uvec ind;

	u1 = Rfast::colMedian(y);
	z = y.each_row() - u1;
	ww = 1 / sqrt(sum(square(z), 1));
	ind = find_finite(ww);
	y = y.rows(ind);
	z = z.rows(ind);
	ww = ww(ind);
	u2 = (ww / accu(ww)) * y;
	while (sum(abs(u2 - u1)) > tol)
	{
		z = y.each_row() - u2;
		u1 = u2;
		ww = 1 / sqrt(sum(square(z), 1));
		if (is_finite(ww.max()))
		{
			u2 = (ww / accu(ww)) * y;
		}
	}
	return u2;
}

RcppExport SEXP Rfast_spat_med(SEXP xSEXP, SEXP tolSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const double>::type tol(tolSEXP);
	__result = spat_med(x, tol);
	return __result;
	END_RCPP
}