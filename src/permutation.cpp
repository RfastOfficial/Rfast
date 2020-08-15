//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

using std::sort;
using std::next_permutation;
using std::prev_permutation;

NumericMatrix permutation(NumericVector X,const unsigned int nperm){
	unsigned int i=0;
	const int n=X.size();
	NumericMatrix F(nperm,n);
	mat FF(F.begin(),nperm,n,false);
	rowvec x(X.begin(),n,false);
	sort(x.begin(),x.end());
	do{
		FF.row(i++)=x;
	} while (next_permutation(x.begin(),x.end()) && i<nperm);
	return F;
}

NumericMatrix permutation_next(NumericVector X,const unsigned int nperm){
  unsigned int i=0;
	const int n=X.size();
	NumericMatrix F(nperm,n);
	mat FF(F.begin(),nperm,n,false);
	rowvec x(X.begin(),n,false);
	do{
		FF.row(i++)=x;
	} while (next_permutation(x.begin(),x.end()) && i<nperm);
    return F(Range(0,i-1),Range(0,n-1));
}

NumericMatrix permutation_prev(NumericVector X,const unsigned int nperm){
  unsigned int i=0;
	const int n=X.size();
	NumericMatrix F(nperm,n);
	mat FF(F.begin(),nperm,n,false);
	rowvec x(X.begin(),n,false);
	do{
		FF.row(i++)=x;
	} while (prev_permutation(x.begin(),x.end()) && i<nperm);
    return F(Range(0,i-1),Range(0,n-1));
}

RcppExport SEXP Rfast_permutation_prev(SEXP xSEXP,SEXP npermSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< const int >::type nperm(npermSEXP);
    __result = permutation_prev(x,nperm);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_permutation_next(SEXP xSEXP,SEXP npermSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< const int >::type nperm(npermSEXP);
    __result = permutation_next(x,nperm);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_permutation(SEXP xSEXP,SEXP npermSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< const int >::type nperm(npermSEXP);
    __result = permutation(x,nperm);
    return __result;
END_RCPP
}