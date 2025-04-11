
//Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

using std::sort;

NumericMatrix rbing(const int n, NumericVector Lam) {
  // lam are the q - 1 non zero eigenvalues
  NumericVector llam_full=clone(Lam);
  sort(llam_full.begin(),llam_full.end(),std::greater<double>()); // sort the eigenvalues in desceding order
  llam_full.push_back(0);
  const int qa = llam_full.size();
  rowvec lam_full(llam_full.begin(),qa,false);
  NumericMatrix x(n,qa);
  mat X(x.begin(),n,qa,false);
  rowvec sigacginv = 1 + 2 * lam_full,SigACG = sqrt( 1 / sigacginv ),y2(qa),yp(qa),y(qa);
  const double qa2 = qa/2.0,tmp = -qa2 * log(qa) + 0.5 * (qa - 1);
  double lratio;
  for (int i=0;i < n;) {
    for(int j=0;j<qa;++j)
      yp[j]=R::rnorm(0,SigACG[j]); 
    y = yp / sqrt( sum(square(yp)) );
    y2 = square(y);
    lratio =  - dot( y2,lam_full ) + tmp + qa2 * log( dot( y2,sigacginv ) );
    if ( log(R::runif(0,1)) < lratio) {
      X.row(i++) = y;
    }
  }
  return x;
}

RcppExport SEXP Rfast_rbing(SEXP nSEXP,SEXP lamSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type n(nSEXP);
    traits::input_parameter< NumericVector >::type lam(lamSEXP);
    __result = rbing(n,lam);
    return __result;
END_RCPP
}