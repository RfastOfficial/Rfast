//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include <cmath>
#include "reg_lib.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


//[[Rcpp::export]]
List colrint_mle(NumericMatrix X, IntegerVector id, const bool ranef, const double tol, const int maxiters, const bool parallel){
  // assume that the function is called with var != 0
  int n = X.nrow(),D = X.ncol(),idmx = max(id);
  mat x0(X.begin(),n,D,false);

  vec ni=Tabulate<vec,IntegerVector>(id,idmx);

  List l;
  mat info(D,3);
  mat ranefmat;
  if(ranef)
    ranefmat = mat(D,idmx);

  if(parallel){
      #ifdef _OPENMP
      #pragma omp parallel
      {
      #endif
        vec x,sx,mx,com,d(2),xminb1,mxminb1,hi2,ni2,oneplnid;
        rowvec tmprowvec;
        double sxy=0,S=0,b1=0,down=0,b2=0,sigma=0;
        mat tmpmat;
        int i=0;
        #ifdef _OPENMP
        #pragma omp for
        #endif
        for(int j = 0; j < D; j++){
          x = x0.col(j);

          sx = group_sum_helper<vec,vec,IntegerVector>(x, id, nullptr,&idmx);
          sxy = sum(x);
          mx = sx/ni;
          com = ni % sx;
          b1 = sxy/n;
          xminb1 = x - b1;
          S = sum(xminb1%xminb1);
          mxminb1 = mx-b1;
          hi2 = mxminb1%mxminb1;
          ni2 = ni%ni;
          d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
          oneplnid = 1 + ni * d[0];
          down = n - d[0] * sum(ni2/(oneplnid));
          b2 = (sxy - d[0] * sum(com/(oneplnid)))/down;
          i = 2;
          while (i++ < maxiters && abs(b2 - b1) > tol) {
            b1 = b2;
            xminb1 = x - b1;
            S = sum(xminb1%xminb1);
            mxminb1 = mx-b1;
            hi2 = mxminb1%mxminb1;
            d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
            oneplnid = 1 + ni * d[0];
            down = n - d[0] * sum(ni2/(oneplnid));
            b2 = (sxy - d[0] * sum(com/(oneplnid)))/down;
          }

          sigma = S/n;
          info(j,1) = sigma/(1 + d[0]);
          info(j,0) = sigma - info(j,1);
          info(j,2) = -0.5 * (d[1] + n * (1.83787706640935-log(n) + 1));
          if (ranef){
            ranefmat.row(j) = conv_to<rowvec>::from((mx - b2) % (d[0] * ni/(oneplnid)));
          }

        }
      #ifdef _OPENMP
      }
      #endif
  }
  else {
    vec x,sx,mx,com,d(2),xminb1,mxminb1,hi2,ni2,oneplnid;
    rowvec tmprowvec;
    double sxy=0,S=0,b1=0,down=0,b2=0,sigma=0;
    mat tmpmat;
    int i=0;

    for(int j = 0; j < D; j++){
      x = x0.col(j);

      sx = group_sum_helper<vec,vec,IntegerVector>(x, id, nullptr,&idmx);
      sxy = sum(x);
      mx = sx/ni;
      com = ni % sx;
      b1 = sxy/n;
      xminb1 = x - b1;
      S = sum(xminb1%xminb1);
      mxminb1 = mx-b1;
      hi2 = mxminb1%mxminb1;
      ni2 = ni%ni;
      d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
      oneplnid = 1 + ni * d[0];
      down = n - d[0] * sum(ni2/(oneplnid));
      b2 = (sxy - d[0] * sum(com/(oneplnid)))/down;
      i = 2;
      while (i++ < maxiters && abs(b2 - b1) > tol) {
        b1 = b2;
        xminb1 = x - b1;
        S = sum(xminb1%xminb1);
        mxminb1 = mx-b1;
        hi2 = mxminb1%mxminb1;
        d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
        oneplnid = 1 + ni * d[0];
        down = n - d[0] * sum(ni2/(oneplnid));
        b2 = (sxy - d[0] * sum(com/(oneplnid)))/down;
      }

      sigma = S/n;
      info(j,1) = sigma/(1 + d[0]);
      info(j,0) = sigma - info(j,1);
      info(j,2) = -0.5 * (d[1] + n * (1.83787706640935-log(n) + 1));
      if (ranef){
        ranefmat.row(j) = conv_to<rowvec>::from((mx - b2) % (d[0] * ni/(oneplnid)));
      }

    }

  }

  l["info"] = info;
  if (ranef){
    l["ranef"] = ranefmat;
  }
  return l;
}

RcppExport SEXP Rfast_colrint_mle(SEXP XSEXP,SEXP idSEXP,SEXP ranefSEXP,SEXP tolSEXP,SEXP maxitersSEXP,SEXP parallelSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< IntegerVector >::type id(idSEXP);
  traits::input_parameter< const bool >::type ranef(ranefSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  __result = colrint_mle(X,id,ranef,tol,maxiters,parallel);
  return __result;
  END_RCPP
}

