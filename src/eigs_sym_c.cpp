#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
List eigs_sym_c(NumericMatrix X,const int k, const bool vectors){
  List l;
  mat x(X.begin(),X.nrow(),X.ncol(),false);
  vec eigval;
  mat eigvec;

  eigs_sym( eigval, eigvec, conv_to<sp_mat>::from(x), k);

  l["values"] = flipud(eigval);
  
  if(vectors){
	l["vectors"] = fliplr(eigvec);
  }
  
  return l;
}


RcppExport SEXP Rfast_eigs_sym_c(SEXP XSEXP,SEXP kSEXP,SEXP vectorsSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type X(XSEXP);
    traits::input_parameter< const int >::type k(kSEXP);
	traits::input_parameter< const bool >::type vectors(vectorsSEXP);
    __result = eigs_sym_c(X,k,vectors);
    return __result;
END_RCPP
}