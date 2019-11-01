//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;

double Norm(NumericMatrix x, const char type) {
	double s=0.0;
	if(type=='F')
		s=std::sqrt(sum_with< square2<double>,NumericMatrix>(x));
	else{
		mat xx(x.begin(),x.nrow(),x.ncol(),false);
		switch(type){
			case 'C':{
				rowvec a=sum(abs(xx),0);
				s=a[a.index_max()];
				break;
			}
			case 'R':{
				colvec a=sum(abs(xx),1);
				s=a[a.index_max()];
				break;
			}
			case 'M':{
				s=xx[xx.index_max()];
				break;
			}
			default:{
				stop("Wrong type. You have to give one of <F,C,R,M>.\n");
			}
		}
	}
	return s;
}


RcppExport SEXP Rfast_Norm(SEXP xSEXP,SEXP typeSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const char  >::type type(typeSEXP);
    __result = wrap(Norm(x,type));
    return __result;
END_RCPP
}