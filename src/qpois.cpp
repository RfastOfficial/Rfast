#include "calc_qpois_regs.h"
#include "mn.h"
#include <RcppArmadillo.h>

arma::vec qpois_regs(arma::mat x, arma::vec y, const double tol, const double ylogy, const double my) {
	return calc_qpois_regs(x, y, tol, ylogy, my);
}


RcppExport SEXP Rfast_qpois_regs(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP,SEXP ylogySEXP,SEXP mySEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< mat >::type x(xSEXP);
    traits::input_parameter< vec >::type y(ySEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const double >::type ylogy(ylogySEXP);
    traits::input_parameter< const double >::type my(mySEXP);
    __result = wrap(qpois_regs(x,y,tol,ylogy,my));
    return __result;
END_RCPP
}

//[[Rcpp::export]]
List qpois_reg(NumericMatrix X, NumericVector Y,const double ylogy,const double tol,const int maxiters){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=pcols;

  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),y(Y.begin(),n,false),m(n),phi(n);
  mat L2,x(X.begin(),n,pcols,false),x_tr(n,pcols);
  double dif;
  b_old(0)=log(mean(y));
  x_tr=x.t();
  int ij =2;
  for(dif=1.0;dif>tol;){
    yhat=x*b_old;
    m=(exp(yhat));
    phi=y-m;
    L1=x_tr*phi;
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
    if(++ij==maxiters)
      break;
  }
  List l;

  l["deviance"]=2.0*(ylogy-sum(y%yhat));
  l["be"]=b_new;
  l["L2"]=L2;
  l["phi"]=sum(arma::square(phi)/m)/(n-pcols);
  l["glm"] = (b_new[d-1]*b_new[d-1])/((sum(arma::square(phi)/m)/(n-pcols))*((mat)solve(L2,mat(d,d,fill::eye),solve_opts::fast))(d-1,d-1));

  return l;
}

RcppExport SEXP Rfast_qpois_reg(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP,SEXP maxitersSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< const double >::type ylogy(ylogySEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    __result = wrap(qpois_reg(x,y,ylogy,tol,maxiters));
    return __result;
END_RCPP
}