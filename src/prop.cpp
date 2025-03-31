#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "mn.h"
#include <cmath>
#include "reg_lib.h"
#include "Rfast.h"

using namespace arma;
using namespace std;

//[[Rcpp::export]]
NumericVector prop_regs(NumericMatrix X, NumericVector Y,const double tol,const std::string varb, const int maxiters){
  const unsigned int n=X.nrow(),pcols=X.ncol();

  unsigned int i=0;
  int ij=0;
  colvec y(Y.begin(),n,false), be(2),expyhat(n),W(n,fill::zeros),x_col(n),x2_col(n),de(n), yhat(n), p(n);
  mat x(X.begin(),n,pcols,false);
  NumericVector F(pcols);
  double d1=0,d2=0,t=0,dera2=0.0,sp=0.0,derb=0.0,dera=0.0,derab=0.0,derb2=0.0;
  const double my=mean(y),sy=my*n,lgmy=log(my/(1-my)),d0 = (sy*log(my)+(n-sy)*log(1-my));//d0 = -2*(sy*log(my)+(n-sy)*log(1-my));

  double W0=my*(1-my);
  double dera20=n*W0,vb;
  colvec de0(n);
  de0=y-my;
  double dera0=0;

  for(i=0;i<pcols;++i){
    d1=d0;
    be[0]=lgmy;
    be[1]=0;
    x_col=x.col(i);
    x2_col=arma::square(x_col);
    derb=sum(de0%x_col);
    derab=W0*sum(x_col);
    derb2=W0*sum(x2_col);
    t=dera20 * derb2 - derab*derab;
    be[0]=be[0]+(derb2 * dera0 - derab * derb)/t;
    be[1]=be[1]+( - derab * dera0 + dera20 * derb )/t;

    yhat = be[0]+be[1]*x_col;
    expyhat=exp(-yhat);
    p = 1 / ( 1 + expyhat );

    d2 = calcDevRes(p,y,yhat);
    ij=2;

    while(ij++<maxiters && (d2-d1)>tol){
      d1=d2;
      W=p%(1-p);
      dera2=sum(W);
      sp=sum(p);

      de=y-p;
      dera=sy-sp;
      derb=sum(de%x_col);
      derab=sum(W%x_col);
      derb2=sum(W%x2_col);
      t=dera2 * derb2 - derab*derab;
      be[0]=be[0]+(derb2 * dera - derab * derb)/t;
      be[1]=be[1]+( - derab * dera + dera2 * derb )/t;
      yhat = be[0]+be[1]*x_col;
      expyhat=exp(-yhat);
      p = 1 / ( 1 + expyhat );

      d2 = calcDevRes(p,y,yhat);
    }

    colvec u2=arma::square(de);

    if(varb=="quasi"){
      double b12=sum(u2%x_col);

      vb=-derab*(-derab*sum(u2)+dera2*b12)/((dera2*derb2-derab*derab)*(dera2*derb2-derab*derab))+dera2*(-derab*b12+sum(u2%x2_col)*dera2)/((dera2*derb2-derab*derab)*(dera2*derb2-derab*derab));
    }else if(varb=="glm"){
      vb=sum(u2/W)/(n-2);
      vb=vb*dera2/(dera2*derb2-derab*derab);
    }else{
      stop("Unsupported varb. Enter \"glm\" or \"quasi\".");
    }
    F[i]=be[1]*be[1]/vb;
  }
  return F;
}



RcppExport SEXP Rfast_prop_regs(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP,SEXP varbSEXP,SEXP maxitersSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const string >::type varb(varbSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    __result = prop_regs(x,y,tol,varb,maxiters);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////////////////////////////////////////


//[[Rcpp::export]]
List prop_reg(NumericMatrix X, NumericVector Y,const double tol,const int maxiters){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=pcols;
  colvec be(d,fill::zeros),yhat(n),expyhat,y(Y.begin(),n,false),W(n,fill::zeros),p(n);
  mat x(X.begin(),n,pcols,false);
  double my = accu(y)/n,d1=(n*my*log(my)+(n-n*my)*log(1-my)),d2;
  be(0)=log(my)-log(1-my);
  mat der=cross_x_y<mat,mat,colvec>(x,y-my);
  mat der2=cross_x_y<mat,mat,colvec>(x,x*my*(1-my));
  be=be+solve(der2,der);
  yhat = x*be;
  expyhat=exp(-yhat);
  p = 1 / (1 + expyhat);
  d2=calcDevRes(p,y,expyhat);
  int i=2;
  for(;d2-d1>tol && i<maxiters;++i){
    d1=d2;
    der=cross_x_y<mat,mat,colvec>(x,y-p);
    W=p%(1-p);
    der2=cross_x_y<mat,mat,colvec>(x,x.each_col()%W);
    be=be+solve(der2,der);
    yhat = x*be;
    expyhat=exp(-yhat);
    p = 1 / (1 + expyhat);
    d2=calcDevRes(p,y,expyhat);
  }
  List l;
  l["p"]=p;
  l["be"]=be;
  l["der2"]=der2;
  l["i"]=i;
  return l;
}

RcppExport SEXP Rfast_prop_reg(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP,SEXP maxitersSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    __result = prop_reg(x,y,tol,maxiters);
    return __result;
END_RCPP
}