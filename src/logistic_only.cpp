//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "mn.h"
#include <cmath>
#include "reg_lib.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

double calc_neg_ll(double *wx, double *expwx, double *y, const int size){
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
NumericVector logistic_only(NumericMatrix X, NumericVector Y,const double tol){
  int maxiters = 100;
  const unsigned int N=X.nrow(), P=X.ncol();

  vec y(Y.begin(),N,false);
  mat x(X.begin(),N,P,false);
  NumericVector F(P);

  double alpha = 1e-4, beta = 0.5, lltol = 1e-06, ttol = 1e-09;
  mat eye(2,2,fill::eye);

  #ifdef _OPENMP
  #pragma omp parallel
  {
  #endif
    vec B(2), nextB, expwxinv = vec(N), pp, u;
    double prevNegLL, negLL,ta;
    mat g, lambda, wx, expwx, der2,tmpX(N,2);
    tmpX.col(0).fill(1);
    int iters;

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(unsigned int i=0;i<P;i++){
      negLL = 0.6931472*N;
      tmpX.col(1) = x.col(i);
      expwxinv.fill(0.5);
      B[0] = 0;
      B[1] = 0;

      for(iters=0; iters<maxiters; iters++){
        g = cross_x_y<mat,mat,vec>(tmpX, expwxinv-y);

        pp = expwxinv % (1 - expwxinv);
        der2 = cross_x_y<mat,mat,vec>(tmpX, tmpX.each_col() % pp);
        u =  solve( der2, g, solve_opts::fast);

        lambda = cross_x_y<mat,mat,vec>(g,u);
        ta = 1/beta;
        // Backtracking line search
        prevNegLL = negLL;
        do{
          ta = ta * beta;
          nextB = B + ta * u;
          wx = tmpX * nextB;
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

      F(i) = 2 * negLL;
    }
  #ifdef _OPENMP
  }
  #endif

  return F;
}

// logistic
RcppExport SEXP Rfast_logistic_only(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type x(xSEXP);
  traits::input_parameter< NumericVector >::type y(ySEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  __result = logistic_only(x,y,tol);
  return __result;
  END_RCPP
}

//[[Rcpp::export]]
NumericMatrix logistic_only_b(NumericMatrix X, NumericVector Y,const double tol){
  int maxiters = 100;
  const unsigned int N=X.nrow(), P=X.ncol();

  vec y(Y.begin(),N,false);
  mat x(X.begin(),N,P,false);
  NumericMatrix F(3,P);

  double alpha = 1e-4, beta = 0.5, lltol = 1e-06, ttol = 1e-09;
  mat eye(2,2,fill::eye);

  #ifdef _OPENMP
  #pragma omp parallel
  {
  #endif
    vec B(2), nextB, expwxinv = vec(N), pp, u;
    double prevNegLL, negLL,ta;
    mat g, lambda, wx, expwx, der2,tmpX(N,2);
    tmpX.col(0).fill(1);
    int iters;

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(unsigned int i=0;i<P;i++){
      negLL = 0.6931472*N;
      tmpX.col(1) = x.col(i);
      expwxinv.fill(0.5);
      B[0] = 0;
      B[1] = 0;

      for(iters=0; iters<maxiters; iters++){
        g = cross_x_y<mat,mat,vec>(tmpX, expwxinv-y);

        pp = expwxinv % (1 - expwxinv);
        der2 = cross_x_y<mat,mat,vec>(tmpX, tmpX.each_col() % pp);
        u =  solve( der2, g, solve_opts::fast);

        lambda = cross_x_y<mat,mat,vec>(g,u);
        ta = 1/beta;
        // Backtracking line search
        prevNegLL = negLL;
        do{
          ta = ta * beta;
          nextB = B + ta * u;
          wx = tmpX * nextB;
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

      F(0,i) = 2 * negLL;
      F(1,i) = -B[0];
      F(2,i) = -B[1];
    }
  #ifdef _OPENMP
  }
  #endif

  return F;
}

RcppExport SEXP Rfast_logistic_only_b(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type x(xSEXP);
  traits::input_parameter< NumericVector >::type y(ySEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  __result = logistic_only_b(x,y,tol);
  return __result;
  END_RCPP
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

NumericVector poisson_only(NumericMatrix X, NumericVector Y,const double ylogy,const double tol){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=2;
  unsigned int i;
  colvec b_old(d),b_new(d),L1(d),yhat(n),y(Y.begin(),n,false);
  mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones),x(X.begin(),n,pcols,false);
  vec m(n),z_col_1(n);
  NumericVector F(pcols);
  double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(mean(y));
  for(i=0;i<pcols;++i){
    b_old[0]=lgmeany;
    b_old[1]=0;
    z_col_1=x.col(i);
    z.col(1)=z_col_1;
    z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
    for(dif=1.0;dif>tol;){
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
    F[i]=2.0*(ylogy-sum(y%yhat));
  }
  return F;
}

// poisson
RcppExport SEXP Rfast_poisson_only(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< const double >::type ylogy(ylogySEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    __result = poisson_only(x,y,ylogy,tol);

    return __result;
END_RCPP
}


//[[Rcpp::export]]
NumericMatrix poisson_only_b(NumericMatrix X, NumericVector Y,double ylogy,const double tol){
  const unsigned int n=X.nrow(),pcols=X.ncol();
  unsigned int i,d=2;
  colvec b_old(d),b_new(d),L1(d),yhat(n),y(Y.begin(),n,false);
  mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n),x(X.begin(),n,pcols,false);
  vec m(n),z_col_1(n);
  NumericMatrix F(3,pcols);
  double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(mean(y));
  for(i=0;i<pcols;++i){
  	b_old(0)=lgmeany;
  	b_old(1)=0;
    z.col(1)=x.col(i);
    z_col_1=z.col(1);
    z_tr=z.t();
    for(dif=1.0;dif>tol;){
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
    F(0,i)=2.0*(ylogy-dot(y,yhat));
    F(1,i)=b_new(0);
    F(2,i)=b_new(1);
  }
  return F;
}

// poisson
RcppExport SEXP Rfast_poisson_only_b(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< double >::type ylogy(ylogySEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    __result = poisson_only_b(x,y,ylogy,tol);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////////////////////////////////////////



NumericMatrix quasi_poisson_only(NumericMatrix X, NumericVector Y, const double ylogy, const double tol,const int maxiters){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=2;
  unsigned int i;
  int ij;

  colvec b_old(d),b_new(d),L1(d),yhat(n);
  vec y(Y.begin(),n,false);
  mat x(X.begin(),n,pcols,false), z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones);
  vec m(n),z_col_1(n);
  NumericMatrix F(2,pcols);
  double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(mean(y));
  for(i=0;i<pcols;++i){
    b_old(0)=lgmeany;
    b_old(1)=0;
    z_col_1=x.col(i);
    z.col(1)=z_col_1;
    z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
    ij=2;
    for(dif=1.0;dif>tol;){
      sm=szm=sz2m=0.0;
      yhat=z*b_old;
      m=(exp(yhat));
      L1=z_tr*(y-m);
      sm=sum(m);
      szm=sum(m%z_col_1);
      sz2m=sum(m%arma::square(z_col_1));
      t=1.0/(sm*sz2m-szm*szm);
      inv_L2.at(0,0)=sz2m*t;
      inv_L2.at(0,1)=inv_L2.at(1,0)=-szm*t;
      inv_L2.at(1,1)=sm*t;
      b_new=b_old+inv_L2*L1;
      dif=sum(abs(b_new-b_old));
      b_old=b_new;
      if(++ij==maxiters)
        break;
    }
    F(0,i)= 2.0*(ylogy-sum(y%yhat));
    F(1,i) = sum(arma::square(y-m)/m)/(n-d);
  }
  return F;
}

RcppExport SEXP Rfast_quasi_poisson_only(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type x(xSEXP);
  traits::input_parameter< NumericVector >::type y(ySEXP);
  traits::input_parameter< const double >::type ylogy(ylogySEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = quasi_poisson_only(x,y,ylogy,tol,maxiters);
  return __result;
  END_RCPP
}
