// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rinternals.h>
#include <R.h>

using namespace Rcpp;

using std::abs;

static void init(double *start,double *end,int s[]){
  *s=s[1]=s[2]=s[3]=0;
  for(;start!=end;++start)
    ++s[(int)abs(*start)];
}

//[[Rcpp::export]]
SEXP odds_helper(SEXP x){
  const int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F=Rf_allocMatrix(INTSXP,4,ncol);
  double *xx=REAL(x),*end=xx+ncol*nrow;
  int *f=INTEGER(F);
  for(;xx!=end;xx+=nrow,f+=4)
    init(xx,xx+nrow,f);
  return F;
}


RcppExport SEXP Rfast_odds_helper(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = odds_helper(x);
    return __result;
END_RCPP
}