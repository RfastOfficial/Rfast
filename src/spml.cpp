///Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include <cmath>
#include "mn.h"
#include "Rfast.h"
#include "reg_lib.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

NumericMatrix spml_regs(NumericMatrix Y, NumericMatrix X0, const double tol, const bool logged, const int maxiters, const bool parallel){
  int n = X0.nrow(),D = X0.ncol(), yD = Y.ncol();
  mat x(X0.begin(),n,D,false), u(n,2), y(Y.begin(), n, yD,false);

  if (yD > 1) {
    u = y;
  } else {
    u.col(0) = cos(y);
    u.col(1) = sin(y);
  }

  vec ci2 = u.col(0)%u.col(0), cisi = u.col(0)%u.col(1), si2 = u.col(1)%u.col(1);

  double f = -0.5, con = 2.506628274631;

  vec one(n,fill::ones);
  double ini = spml_mle2(u,ci2,cisi,si2,n,1e-09,maxiters);
  NumericMatrix ret(D,2);

  if(parallel){
      #ifdef _OPENMP
     #pragma omp parallel
     {
     #endif
    double sumX, lik1, lik2;
    int i;
    mat xx(2,2,fill::zeros),X(n,2), B,mu,a11,a12,a22,der2(4,4),XX,der(4,1);
    vec tau,ptau,psit,rat,psit2,slv;
    mat::iterator tmpit = der2.begin(),it;
    X.col(0) = one;

    xx(0) = n;
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(int j = 0; j < D; j++){
      X.col(1) = x.col(j);
      sumX = sum(X.col(1));
      xx(1) = sumX;
      xx(2) = sumX;
      xx(3) = sum(X.col(1)%X.col(1));
      XX = solve(xx, X.t());
      B =  XX * u;
      mu = X * B;
      tau = sum(u % mu,1);
      ptau = pnormc(tau);

      lik1 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n);

      rat = ptau / ( exp(f * (tau%tau))/con + tau%ptau );
      psit = tau + rat;
      psit2 = 2 -((tau+rat)%rat);
      a11 = cross_x_y<mat,mat,vec>(X, u.each_col()%psit-mu);
      der[0] = a11[0],der[1] = a11[1],der[2] = a11[2],der[3] = a11[3];
      a11 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%ci2 - 1));
      a12 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%cisi));
      a22 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%si2 - 1));
      it = tmpit;

      // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
      (*(it)) = a11[0], (*(++it)) = a11[2], (*(++it)) = a12[0], (*(++it)) = a12[2];
      (*(++it)) = a11[1], (*(++it)) = a11[3], (*(++it)) = a12[1], (*(++it)) = a12[3];
      (*(++it)) = a12[0], (*(++it)) = a12[2], (*(++it)) = a22[0], (*(++it)) = a22[2];
      (*(++it)) = a12[1], (*(++it)) = a12[3], (*(++it)) = a22[1], (*(++it)) = a22[3];

      slv = solve(der2,der);
      B[0] = B[0]-slv[0];
      B[1] = B[1]-slv[1];
      B[2] = B[2]-slv[2];
      B[3] = B[3]-slv[3];

      mu = X*B;

      tau = sum(u % mu,1);
      ptau = pnormc(tau);

      lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n);

      i=2;
      while(i++<maxiters && (lik2 - lik1) > tol){
        lik1 = lik2;
        rat = ptau / ( exp(f * (tau%tau))/con + tau%ptau );
        psit = tau + rat;
        psit2 = 2 - tau % rat - rat%rat;
        a11 = cross_x_y<mat,mat,vec>(X, u.each_col()%psit-mu);
        der[0] = a11[0],der[1] = a11[1],der[2] = a11[2],der[3] = a11[3];
        a11 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%ci2 - 1));
        a12 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%cisi));
        a22 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%si2 - 1));
        it = tmpit;

        // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
        (*(it)) = a11[0], (*(++it)) = a11[2], (*(++it)) = a12[0], (*(++it)) = a12[2];
        (*(++it)) = a11[1], (*(++it)) = a11[3], (*(++it)) = a12[1], (*(++it)) = a12[3];
        (*(++it)) = a12[0], (*(++it)) = a12[2], (*(++it)) = a22[0], (*(++it)) = a22[2];
        (*(++it)) = a12[1], (*(++it)) = a12[3], (*(++it)) = a22[1], (*(++it)) = a22[3];

        slv = solve(der2,der);
        B[0] = B[0]-slv[0];
        B[1] = B[1]-slv[1];
        B[2] = B[2]-slv[2];
        B[3] = B[3]-slv[3];

        mu = X*B;

        tau = sum(u % mu,1);
        ptau = pnormc(tau);

        lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n);
      }

      ret(j,0) = 2 * (lik2 - n * 1.83787706640935 - ini);
      ret(j,1) = R::pchisq(ret(j,0), 2, false, logged);
    }
    #ifdef _OPENMP
    }
    #endif
  }
  else{
    double sumX, lik1, lik2;
    int i;
    mat xx(2,2,fill::zeros),X(n,2), B,mu,a11,a12,a22,der2(4,4),XX,der(4,1);
    vec tau,ptau,psit,rat,psit2,slv;
    mat::iterator tmpit = der2.begin(),it;
    X.col(0) = one;

    xx(0) = n;

    for(int j = 0; j < D; j++){
      X.col(1) = x.col(j);
      sumX = sum(X.col(1));
      xx(1) = sumX;
      xx(2) = sumX;
      xx(3) = sum(X.col(1)%X.col(1));
      XX = solve(xx, X.t());
      B =  XX * u;
      mu = X * B;
      tau = sum(u % mu,1);
      ptau = pnormc(tau);

      lik1 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n);

      rat = ptau / ( exp(f * (tau%tau))/con + tau%ptau );
      psit = tau + rat;
      psit2 = 2 - tau % rat - rat%rat;
      a11 = cross_x_y<mat,mat,vec>(X, u.each_col()%psit-mu);
      der[0] = a11[0],der[1] = a11[1],der[2] = a11[2],der[3] = a11[3];
      a11 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%ci2 - 1));
      a12 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%cisi));
      a22 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%si2 - 1));
      it = tmpit;

      // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
      (*(it)) = a11[0], (*(++it)) = a11[2], (*(++it)) = a12[0], (*(++it)) = a12[2];
      (*(++it)) = a11[1], (*(++it)) = a11[3], (*(++it)) = a12[1], (*(++it)) = a12[3];
      (*(++it)) = a12[0], (*(++it)) = a12[2], (*(++it)) = a22[0], (*(++it)) = a22[2];
      (*(++it)) = a12[1], (*(++it)) = a12[3], (*(++it)) = a22[1], (*(++it)) = a22[3];

      slv = solve(der2,der);
      B[0] = B[0]-slv[0];
      B[1] = B[1]-slv[1];
      B[2] = B[2]-slv[2];
      B[3] = B[3]-slv[3];

      mu = X*B;

      tau = sum(u % mu,1);
      ptau = pnormc(tau);

      lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n);

      i=2;
      while(i++<maxiters && (lik2 - lik1) > tol){
        lik1 = lik2;
        rat = ptau / ( exp(f * (tau%tau))/con + tau%ptau );
        psit = tau + rat;
        psit2 = 2 - tau % rat - rat%rat;
        a11 = cross_x_y<mat,mat,vec>(X, u.each_col()%psit-mu);
        der[0] = a11[0],der[1] = a11[1],der[2] = a11[2],der[3] = a11[3];
        a11 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%ci2 - 1));
        a12 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%cisi));
        a22 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%si2 - 1));
        it = tmpit;

        // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
        (*(it)) = a11[0], (*(++it)) = a11[2], (*(++it)) = a12[0], (*(++it)) = a12[2];
        (*(++it)) = a11[1], (*(++it)) = a11[3], (*(++it)) = a12[1], (*(++it)) = a12[3];
        (*(++it)) = a12[0], (*(++it)) = a12[2], (*(++it)) = a22[0], (*(++it)) = a22[2];
        (*(++it)) = a12[1], (*(++it)) = a12[3], (*(++it)) = a22[1], (*(++it)) = a22[3];

        slv = solve(der2,der);
        B[0] = B[0]-slv[0];
        B[1] = B[1]-slv[1];
        B[2] = B[2]-slv[2];
        B[3] = B[3]-slv[3];

        mu = X*B;

        tau = sum(u % mu,1);
        ptau = pnormc(tau);

        lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n);
      }
      ret(j,0) = 2 * (lik2 - n * 1.83787706640935 - ini);
      ret(j,1) = R::pchisq(ret(j,0), 2, false, logged);
    }
  }

  return ret;
}

RcppExport SEXP Rfast_spml_regs(SEXP YSEXP,SEXP X0SEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP maxitersSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    traits::input_parameter< NumericMatrix >::type X0(X0SEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const bool >::type logged(loggedSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    traits::input_parameter< const int >::type parallel(parallelSEXP);
    __result = spml_regs(Y,X0,tol,logged,maxiters,parallel);
    return __result;
END_RCPP
}



/////////////////////////////////////////////////////////////////////////////////////////////////


//[[Rcpp::export]]
List spml_reg(NumericMatrix Y, NumericMatrix X, const double tol, const bool seb, const int maxiters){
  int n = X.nrow(), D = X.ncol(), yD = Y.ncol();
  List l;
  vec ci, si;
  mat x(X.begin(),n,D,false), y(Y.begin(), n, yD,false), u(n,2), tx = x.t(), XX = solve(tx*x,tx);

  if (yD > 1) {
    u = y;
  } else {
    u.col(0) = cos(y);
    u.col(1) = sin(y);
  }

  mat be = XX*u, mu = x*be;

  vec tau = sum(u%mu,1), ptau = pnormc(tau), ci2 = u.col(0)%u.col(0), cisi = u.col(0)%u.col(1), si2 = u.col(1)%u.col(1);

  double f= -0.5, con = 2.506628274631, lik1 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n), lik2;
  vec rat, psit2, psit = tau + ptau / ( exp(f * tau%tau)/con + tau % ptau );
  be = XX * (u.each_col() % psit);
  mu = mu.each_col()% psit;
  tau = sum(u%mu,1);
  ptau = pnormc(tau);
  lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n);
  mat der2, a11,a12,a22;

  mat der, slv;
  int i = 2;

  while(i++<maxiters && lik2 - lik1 > tol){
    lik1 = lik2;

    rat = ptau / ( exp(f * tau%tau)/con + tau % ptau );
    psit = tau + rat;

    psit2 = 2 -((tau+rat)%rat);
    der = cross_x_y<mat,mat,vec>(x, u.each_col()%psit-mu);
    der.reshape(2*D,1);



    a11 = cross_x_y<mat,mat,vec>(x, x.each_col() % (psit2%ci2 - 1));
    a12 = cross_x_y<mat,mat,vec>(x, x.each_col() % (psit2%cisi));
    a22 = cross_x_y<mat,mat,vec>(x, x.each_col() % (psit2%si2 - 1));

    der2 = join_cols(join_rows(a11,a12),join_rows(a12,a22));

    // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
    /*
    (*(tmpit)) = a11[0], (*(++tmpit)) = a11[2], (*(++tmpit)) = a12[0], (*(++tmpit)) = a12[2];
    (*(++tmpit)) = a11[1], (*(++tmpit)) = a11[3], (*(++tmpit)) = a12[1], (*(++tmpit)) = a12[3];
    (*(++tmpit)) = a12[0], (*(++tmpit)) = a12[2], (*(++tmpit)) = a22[0], (*(++tmpit)) = a22[2];
    (*(++tmpit)) = a12[1], (*(++tmpit)) = a12[3], (*(++tmpit)) = a22[1], (*(++tmpit)) = a22[3];
     */

    slv = solve(der2, der);
    slv.reshape(D,2);
    be = be - slv;

    mu = x * be;
    tau = sum(u%mu,1);
    ptau = pnormc(tau);

    lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),tau.begin(),ptau.begin(),n);
  }

  l["iters"] = i;
  l["loglik"] = lik2 - n*1.83787706641;
  l["be"] = be;
  if ( seb ) {
    mat seb = sqrt(((mat)solve( -der2, eye<mat>(2*D,2*D) )).diag());

    seb.reshape(D,2);
    l["seb"] =  seb;
  }
  return l;
}

RcppExport SEXP Rfast_spml_reg(SEXP YSEXP, SEXP XSEXP,SEXP tolSEXP,SEXP sebSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type Y(YSEXP);
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const bool >::type seb(sebSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = wrap(spml_reg(Y,X,tol,seb,maxiters));
  return __result;
  END_RCPP
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////


List spml_mle(NumericMatrix X, const double tol, const int maxiters){
  int n = X.nrow(), yD = X.ncol();
  mat x(X.begin(),n, yD,false);
  List l;
  mat u(n,2);
  if (yD > 1) {
    u = x;
  } else {
    u.col(0) = cos(x);
    u.col(1) = sin(x);
  }

  vec ci2 = u.col(0)%u.col(0),cisi = u.col(0)%u.col(1), si2 = u.col(1)%u.col(1);

  vec su(2);
  su[0] = sum(u.col(0)), su[1] = sum(u.col(1));

  double nR = sqrt(sum(su%su)), kappa = vmf_mle2(nR, n, tol, maxiters);

  vec mu = su*(1/nR);

  vec mu1 = mu*kappa;
  double f = -0.5, con = 2.506628274631;
  vec tau = u*mu1, ptau = pnormc(tau);

  vec rat = ptau/(exp(f * tau%tau)/con + tau % ptau);

  vec psit = tau + rat;
  vec psit2 = 2 - rat%(tau + rat);

  vec der(2);
  der[0] = sum(u.col(0)%psit) - n * mu1[0];
  der[1] = sum(u.col(1)%psit) - n * mu1[1];

  double dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;

  double down = dera2 * derb2 - derab*derab;
  vec mu2(2);
  mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
  mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

  int i = 2;
  while (i++<maxiters && sum_with<abs, vec>(mu2 - mu1) > tol) {
    mu1 = mu2;
    tau = u*mu1;
    ptau = pnormc(tau);
    rat = ptau/(exp(f * (tau%tau))/con + tau % ptau);
    psit = tau + rat;
    psit2 = 2 - rat%(tau + rat);

    der[0] = sum(u.col(0)%psit) - n * mu1[0];
    der[1] = sum(u.col(1)%psit) - n * mu1[1];
    dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;
    down = dera2 * derb2 - derab*derab;
    mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
    mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

  }

  l["iters"] = i-1;

  double gam = mu2[0]*mu2[0]+mu2[1]*mu2[1];

  l["loglik"] = -n * (0.5  * gam + 1.83787706640935) + sum_with<log1p, colvec>(((tau % ptau) * con)/exp(f*tau%tau));

  l["gamma"] = gam;
  l["mu"] = conv_to<rowvec>::from(mu2);

  return l;
}

RcppExport SEXP Rfast_spml_mle(SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = wrap(spml_mle(X,tol,maxiters));
  return __result;
  END_RCPP
}

////////////////////////////////////////////////////////////////////////////////////////////////

