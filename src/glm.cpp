//Author: Manos Papadakis
// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;

static double calc_neg_ll(double *wx, double *expwx, double *y, const int size){
  double sum = 0.0;
  double *wit = wx, *yit = y;
  for(int i=0;i<size;i++,wit++,yit++){
    if(*wit<=30)
      sum+=(*yit-1)*(*wit)+log(expwx[i]);
    else
      sum+=(*yit)*(*wit);
  }
  return sum;
}

//[[Rcpp::export]]
List glm_logistic(NumericMatrix X, NumericVector Y, const double tol = 1e-06, const int maxiters = 100){
  int N = X.nrow(),P = X.ncol();
  mat x(X.begin(), N, P,false);
  vec y(Y.begin(), N,false);
  List l;


  double ta, alpha = 1e-4, beta = 0.5, lltol = 1e-06, ttol = 1e-09;
  vec B(P, fill::zeros), nextB, expwxinv = vec(N,fill::zeros)+0.5, pp, u;
  double prevNegLL, negLL = 0.6931472*N;

  //0->newton, 1->CONJUGATE, 2->GRADIENT

  int iters;
  mat g, lambda, eye(P,P,fill::eye), wx, expwx, der2;

  for(iters=0; iters<maxiters; iters++){
    g = cross_x_y<mat,mat,vec>(x, expwxinv-y);

    pp = expwxinv % (1 - expwxinv);
    der2 = cross_x_y<mat,mat,vec>(x, x.each_col() % pp);
    u =  solve( der2, g, solve_opts::fast);

    lambda = cross_x_y<mat,mat,vec>(g,u);

    // Backtracking line search
    ta = 1/beta;
    prevNegLL = negLL;

    do{
      ta = ta * beta;
      nextB = B + ta * u;
      wx = x * nextB;
      expwx = 1 + exp(wx);
      negLL = calc_neg_ll(&wx[0], &expwx[0], &y[0], N);
    } while(negLL > prevNegLL + alpha * ta * lambda[0] && ta > ttol);

    B = nextB;
    if ( std::isinf(negLL) || lambda[0]*ta/2 < tol || prevNegLL - negLL < lltol) {
      if ( NumericVector::is_na(negLL)) {
        Rcout<<"Infinity found"<<endl;
      }
      break;
    }
    expwxinv = 1/expwx;
  }
  double dev = 2 * negLL;
  if ( std::isinf(dev) )
    dev = 1e+308;

  l["iter"] = iters+1;
  l["be"] = -B;
  l["deviance"] = dev;
  l["der2"] = der2;
  
  return l;
}

RcppExport SEXP Rfast_glm_logistic(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP,SEXP maxitersSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    __result = glm_logistic(x,y,tol,maxiters);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////////////////


List glm_poisson(NumericMatrix X, NumericVector Y,const double ylogy,const double tol){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=pcols;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),y(Y.begin(),n,false),m(n);
  mat L2,x(X.begin(),n,pcols,false),x_tr(n,pcols);
  double dif;
  b_old(0)=log(mean(y));
  x_tr=x.t();
  for(dif=1.0;dif>tol;){
    yhat=x*b_old;
    m=exp(yhat);
    L1=x_tr*(y-m);
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
  }
  List l;
  l["deviance"]=2.0*(ylogy-sum(y%yhat));
  l["be"]=b_new;
  l["L2"]=L2;
  return l;
}

RcppExport SEXP Rfast_glm_poisson(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< const double >::type ylogy(ylogySEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    __result = glm_poisson(x,y,ylogy,tol);
    return __result;
END_RCPP
}

////////////////////////////////////////////////////////////////////////
