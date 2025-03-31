//Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <cmath>
#include "reg_lib.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
NumericMatrix colweibull_mle(NumericMatrix X, const double tol, const int maxiters,const bool parallel){
  int n = X.nrow(),d = X.ncol();
  mat x0(X.begin(),n,d,false);
  NumericMatrix ret(d,3);

  if(parallel){
      #ifdef _OPENMP
      #pragma omp parallel if(parallel)
      {
      #endif
        vec x,lx,lx2,y;
        double mlx,co,sy,fb,fb2,b1,b2;
        int i;
        #ifdef _OPENMP
        #pragma omp for
        #endif
        for(int j = 0; j<d;j++){
          x = x0.col(j), lx = log(x),lx2 = lx%lx, y = x;
          mlx = sum(lx)/n, co = sum(y%lx),sy = sum(y), fb = 1+mlx-co/sy,fb2 = -1 -(sum(y%lx2)*sy-co*co)/(sy*sy);
          b1 = 1,b2 = 1 - fb/fb2;
          i=2;
          while (++i<maxiters && sum(abs(b2 - b1)) > tol) {
            b1 = b2;
            my_pow2(x,&y[0],b1,n);
            co = sum(y % lx);
            sy = sum(y);
            fb = 1/b1 + mlx - co/sy;
            fb2 = -1/(b1*b1) - (sum(y % lx2) * sy - co*co)/(sy*sy);
            b2 = b1 - fb/fb2;
          }

          ret(j,0) = b2;
          ret(j,1) = pow(sy/n,1/b2);
          my_pow2(conv_to<vec>::from(x/ret(j,1) ),&y[0],b2,n);
          ret(j,2) = n * log(b2) - n * b2 * log(ret(j,1) ) + (b2 - 1) * n * mlx - sum(y);
        }
      #ifdef _OPENMP
      }
      #endif
  }
  else{
    vec x,lx,lx2,y;
    double mlx,co,sy,fb,fb2,b1,b2;
    int i;
    for(int j = 0; j<d;j++){
      x = x0.col(j), lx = log(x),lx2 = lx%lx, y = x;
      mlx = sum(lx)/n, co = sum(y%lx),sy = sum(y), fb = 1+mlx-co/sy,fb2 = -1 -(sum(y%lx2)*sy-co*co)/(sy*sy);
      b1 = 1,b2 = 1 - fb/fb2;
      i=2;
      while (++i<maxiters && sum(abs(b2 - b1)) > tol) {
        b1 = b2;
        my_pow2(x,&y[0],b1,n);
        co = sum(y % lx);
        sy = sum(y);
        fb = 1/b1 + mlx - co/sy;
        fb2 = -1/(b1*b1) - (sum(y % lx2) * sy - co*co)/(sy*sy);
        b2 = b1 - fb/fb2;
      }
      ret(j,0) = b2;
      ret(j,1) = pow(sy/n,1/b2);
      my_pow2(conv_to<vec>::from(x/ret(j,1) ),&y[0],b2,n);
      ret(j,2) = n * log(b2) - n * b2 * log(ret(j,1) ) + (b2 - 1) * n * mlx - sum(y);
    }
  }

  return ret;
}



RcppExport SEXP Rfast_colweibull_mle(SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type X(XSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = colweibull_mle(X,tol,maxiters,parallel);
    return __result;
END_RCPP
}
