//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "Rfast.h"
#include "mn.h"

using namespace Rcpp;
using namespace arma;

NumericVector diri_nr_type2(colvec a1,NumericVector A2,colvec ma,int p,const double tol){
  mat f2(p,p),slv;
  NumericVector AA2=clone(A2);
  double sa;
  colvec f,a2(AA2.begin(),AA2.size(),false);
  while ( sum( abs( a2 - a1 ) ) > tol ) {
    a1 = a2;
    sa = sum(a1);
    f = ma - foreach<digamma,colvec>(a1) + digamma(sa);
    f2.fill(trigamma(sa));
    f2.diag() = f2.diag() - foreach<trigamma,colvec>(a1);
    slv = solve(f2, f);
    a2 = a1 - slv.each_col();
  }
  return AA2;
}

// part of the Compositional::diri.nr
RcppExport SEXP Rfast_diri_nr_type2(SEXP a1SEXP,SEXP a2SEXP,SEXP maSEXP,SEXP pSEXP,SEXP tolSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< colvec >::type a1(a1SEXP);
    traits::input_parameter< NumericVector >::type a2(a2SEXP);
    traits::input_parameter< colvec >::type ma(maSEXP);
    traits::input_parameter< int >::type p(pSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    __result = wrap(diri_nr_type2(a1,a2,ma,p,tol));
    return __result;
END_RCPP
}