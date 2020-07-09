// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double total_dista(NumericMatrix Xnew, NumericMatrix X,const bool sqr) {
	const int n=X.ncol(),nu=Xnew.ncol();
	mat xnew(Xnew.begin(),Xnew.nrow(),nu,false),x(X.begin(),X.nrow(),n,false);
	double a=0.0;
    if(sqr){
    	for(int i=0;i<nu;++i){
        	a+=sum_with< square2<double>,mat >(x.each_col() - xnew.col(i));
    	}
    }else{
        for(int i=0;i<nu;++i){
            a+=sum_with<std::sqrt,rowvec>(sum(square(x.each_col() - xnew.col(i)),0));
        }
    }
  	return a;
}

RcppExport SEXP Rfast_total_dista(SEXP xSEXP,SEXP ySEXP,SEXP sqrSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericMatrix >::type y(ySEXP);
    traits::input_parameter< const bool >::type sqr(sqrSEXP);
    __result = total_dista(x,y,sqr);
    return __result;
END_RCPP
}