#include "mn.h"
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

arma::vec qpois_regs(arma::mat x,arma::vec y,const double tol,const double ylogy,const double my) {
  const unsigned int n=x.n_rows,pcols=x.n_cols,d=2;
  unsigned int i;
  colvec b_old(d),b_new(d),L1(d),yhat(n);
  mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones);
  vec m(n),z_col_1(n);
  arma::colvec f(pcols);
  double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(my);
  for(i=0;i<pcols;++i){
    b_old(0)=lgmeany;
    b_old(1)=0;
    z_col_1=x.col(i);
    z.col(1)=z_col_1;
    z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
    for(dif=1.0;dif > 0.000000001 && dif > tol;){
      sm=szm=sz2m=0.0;
      yhat=z*b_old;
      m=exp(yhat);
      L1=z_tr*(y-m);
      sm=sum(m);
      szm=sum(m%z_col_1);
      sz2m=sum(m%square(z_col_1));
      t=1.0/(sm*sz2m-szm*szm);
      inv_L2.at(0,0)=sz2m*t;
      inv_L2.at(0,1)=inv_L2.at(1,0)=-szm*t;
      inv_L2.at(1,1)=sm*t;
      b_new=b_old+inv_L2*L1;
      dif=sum(abs(b_new-b_old));
      b_old=b_new;
    }
  	const double phi=arma::sum(arma::pow(y-m,2)/m)/(n-d); 
	f[i]=(b_new[1]*b_new[1])/(phi*inv_L2.at(1,1));
  }
  return f;
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
    __result = qpois_regs(x,y,tol,ylogy,my);
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
    __result = qpois_reg(x,y,ylogy,tol,maxiters);
    return __result;
END_RCPP
}