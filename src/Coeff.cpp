
#include <string>
#include "Coeff.h"
#include "mn.h"

using Rcpp::NumericMatrix;
using std::string;
using namespace arma;

namespace Coeff
{
    //[[Rcpp::export]]
    NumericMatrix bhattacharyya(NumericMatrix x)
    {
        const int ncl = x.ncol(), nrw = x.nrow();
        mat xx(x.begin(), nrw, ncl, false);
        NumericMatrix f(ncl, ncl);
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
        colvec xv(nrw);
        double a;
        int i, j;

        for (i = 0; i < ncl - 1; ++i)
        {
            xv = sqrt_xx.col(i);
            for (j = i + 1; j < ncl; ++j)
            {
                a = Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j));
                f(i, j) = a;
                f(j, i) = a;
            }
        }
        return f;
    }
}

//[[Rcpp::export]]
NumericMatrix coeff(NumericMatrix x, const string method)
{
    if (method == "bhattacharyya")
    {
        return Coeff::bhattacharyya(x);
    }
    stop("Unsupported Method: %s", method);
}

RcppExport SEXP Rfast_coeff(SEXP xSEXP, SEXP methodSEXP)
{
    BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter<NumericMatrix>::type x(xSEXP);
    traits::input_parameter<const string>::type method(methodSEXP);
    __result = coeff(x, method);
    return __result;
    END_RCPP
}

namespace CoeffVector
{
    //[[Rcpp::export]]
    NumericVector bhattacharyya(NumericMatrix x)
    {
        const int ncl = x.ncol(), nrw = x.nrow();
        mat xx(x.begin(), nrw, ncl, false);
        NumericVector f(proper_size(nrw, ncl));
		mat sqrt_xx(nrw, ncl, fill::none);
		fill_with<std::sqrt, double *, double *>(xx.begin(), xx.end(), sqrt_xx.begin());
        colvec xv(nrw);
        int i, j, k = 0;
        for (i = 0; i < ncl - 1; ++i)
        {
            xv = sqrt_xx.col(i);
            for (j = i + 1; j < ncl; ++j, ++k)
            {
                f[k] = Coeff::bhattacharyya<false>(xv, sqrt_xx.col(j));
            }
        }
        return f;
    }
}

//[[Rcpp::export]]
NumericVector coeff_vec(NumericMatrix x, const string method)
{
    if (method == "bhattacharyya")
    {
        return CoeffVector::bhattacharyya(x);
    }
    stop("Unsupported Method: %s", method);
}

RcppExport SEXP Rfast_coeff_vec(SEXP xSEXP, SEXP methodSEXP)
{
    BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter<NumericMatrix>::type x(xSEXP);
    traits::input_parameter<const string>::type method(methodSEXP);
    __result = coeff_vec(x, method);
    return __result;
    END_RCPP
}
