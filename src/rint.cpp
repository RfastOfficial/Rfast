//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "Rfast.h"
#include <cmath>
#include "reg_lib.h"
#include "mn.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
List rint_reg(NumericMatrix X, NumericVector Y, IntegerVector id, const double tol, const bool ranef, const int maxiters){
  int n = X.nrow(), p = X.ncol(), idmx,idmn;
  maximum<int>(id.begin(),id.end(),idmx);
  minimum<int>(id.begin(),id.end(),idmn);
  mat x(X.begin(), n,p,false),xx(p,p),sx(idmx,p),sxy(p,1),mx(idmx,p);
  vec y(Y.begin(),n,false),my(idmx);

  double logpitimes2 = 1.83787706640935,logn = log(n);
  vec ni=Tabulate<vec,IntegerVector>(id,idmx);

  xx = cross_x<mat,mat>(x);
  for(int i=0;i<p;i++)
    sx.col(i) = group_sum_helper<vec,vec,IntegerVector>(x.col(i), id, &idmn,&idmx);
  sxy = cross_x_y<mat,mat,vec>(x,y);
  colvec sy = group_sum_helper<colvec,vec,IntegerVector>(y, id, &idmn,&idmx);
  mx = sx.each_col()/ni;
  my = sy/ni;

  mat b1 = solve(xx,sxy,solve_opts::fast);
  vec tmp = y - x*b1;
  double S = sum_with<square2<double>, vec>(tmp);
  vec tmp2 = my-mx*b1;
  vec hi2 = tmp2%tmp2;
  vec ni2 = ni%ni;

  vec d(2);
  d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
  vec oneplnid = 1+ni*d(0);
  vec b2 = solve(xx - d(0)* cross_x_y<mat,mat,vec>(sx.each_col()/oneplnid, sx), sxy -
    d(0) * cross_x_y<mat,mat,vec>(sx, sy/oneplnid),solve_opts::fast);
  int i = 2;

  while(i++<maxiters && sum(abs(b2-b1.col(0))) > tol) {
    b1.col(0) = b2;

    tmp = y - x*b1;
	S = sum_with<square2<double>, vec>(tmp);
    tmp2 = my-mx*b1;
    hi2 = tmp2%tmp2;

    d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
    oneplnid = 1+ni*d(0);
    b2 = solve(xx - d(0) * cross_x_y<mat,mat,vec>(sx.each_col()/oneplnid, sx), sxy -
      d(0) * cross_x_y<mat,mat,vec>(sx, sy/oneplnid),solve_opts::fast);
  }

  List l;
  NumericVector info(6);
  info(0) = i-1;
  info(2) = (S-d(0)*sum(ni2%hi2/oneplnid))/n;
  info(1) = d(0)*info(2);
  info(3) = -0.5 * d(1)-0.5*n*(logpitimes2-logn+1);
  info(4) = -2 * info(3);
  info(5) = info(4) + (p + 2) * logn;
  l["info"] = info;
  l["be"] = b2;
  l["seb"] = sqrt(((mat)solve(xx-d(0)*cross_x_y<mat,mat,vec>(sx.each_col()/oneplnid,sx),eye(p,p),solve_opts::fast)).diag()*info(2));

  if(ranef){
    mat er = y - x * (conv_to<colvec>::from(b2));
    l["ranef"] =  d[0] * ni/(oneplnid) % group_sum_helper<vec,vec,IntegerVector>(er.col(0), id, &idmn,&idmx)/ni;
  }
  return l;
}

RcppExport SEXP Rfast_rint_reg(SEXP XSEXP,SEXP YSEXP,SEXP idSEXP,SEXP tolSEXP,SEXP ranefSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< NumericVector >::type Y(YSEXP);
  traits::input_parameter< IntegerVector >::type id(idSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const bool >::type ranef(ranefSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = rint_reg(X,Y,id,tol,ranef,maxiters);
  return __result;
  END_RCPP
}

/////////////////////////////////////////////////////////////////////////


//[[Rcpp::export]]
NumericMatrix rint_regs(NumericMatrix X, NumericVector Y, IntegerVector id, const double tol,
                        const bool logged, const bool parallel, const int maxiters){
  //if(ret != 1 && ret != 2 and ret != 3)
  //  stop("Invalid return option, ret should be 1, 2 or 3\n");
  //1 - stat,pval, 2- all, 3 bic
  int n = X.nrow(), D = X.ncol(), idmx,idmn;
  mat x(X.begin(), n,D,false);
  vec y(Y.begin(),n,false);

  maximum<int>(id.begin(),id.end(),idmx);
  minimum<int>(id.begin(),id.end(),idmn);

  int ret = 1;

  vec ni=Tabulate<vec,IntegerVector>(id,idmx);

  vec ni2 = ni%ni;
  colvec sy = group_sum_helper<colvec,vec,IntegerVector>(y, id, &idmn,&idmx);
  vec my = sy/ni;
  double Sy = sum(sy),logn = log(n);
  vec r = conv_to<vec>::from(cov(y,x));
  double mesi = Sy/n;
  vec xs = conv_to<vec>::from(sum(x));
  vec xs2 = conv_to<vec>::from(sum(x%x));
  vec vx = (xs2 - (xs%xs)/n)/(n - 1);
  vec b(D);
  b = r/vx;
  vec a(D);
  a = mesi - b % xs/n;
  mat be(D,2);
  be.col(0) = a;
  be.col(1) = b;

  NumericMatrix mymat;

  if(ret==1){
    mymat = NumericMatrix(D,2);
  }
  else if(ret==2){
    mymat = NumericMatrix(D,3);
  }
  else{
    mymat = NumericMatrix(D,1);
  }

  if(parallel){
  #ifdef _OPENMP
  #pragma omp parallel
  {
  #endif
    vec Xi(n), sxy(2), b1(2), tmpvec(n), tmpvec2(idmx), hi2(idmx),b2(2),B(2), mx(idmx);
    mat sx(idmx,2), temptcom(idmx,2), tcom(2,idmx), A(2,2);
    vec oneplnid;
    sx.col(0) = ni;
    sx.col(1) = ni;
    int ij=0;
    sxy[0] = Sy;
    double S=0,down=0, se=0,seb=0,info1=0,info2=0,stat,bic;
    vec d(2);
    mat  xx(2,2);
    xx(0,0) = n;
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(int i = 0; i < D; i++) {
      Xi = x.col(i);
      xx(0,1) = xs[i];
      xx(1,0) = xs[i];
      xx(1,1) = xs2[i];
      sx.col(1) = group_sum_helper<vec,vec,IntegerVector>(Xi, id, &idmn,&idmx);
      sxy[1] = sum(Xi % y);
      mx = sx.col(1)/ni;
      b1[0] = be.row(i)[0];
      b1[1] = be.row(i)[1];
      tmpvec = y - b1(0) - b1(1) * Xi;
      S = sum_with<square2<double>, vec>(tmpvec);
      tmpvec2 = my - b1(0) - b1(1) * mx;
      hi2 = tmpvec2%tmpvec2;
      d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
      oneplnid = (1+ ni * d[0]);
      temptcom.col(0) = sx.col(0)/oneplnid;
      temptcom.col(1) = sx.col(1)/oneplnid;
      tcom = -d[0] * temptcom.t();

      A = xx + tcom * sx;

      B = sxy + tcom * sy;

      down = A(0,0) * A(1,1) - A(0,1)*A(0,1);
      b2(0) = (A(1,1) * B(0) - A(0,1) * B(1))/down;
      b2(1) = (- A(0,1) * B(0) + A(0,0) * B(1))/down;
      ij = 2;
      while(ij++<maxiters && sum_with<std::abs,vec>(b2 - b1)>tol){
        b1 = b2;
        tmpvec = y - b1(0) - b1(1) * Xi;
        S = sum_with< square2<double>, vec>(tmpvec);
        tmpvec2 = my - b1(0) - b1(1) * mx;
        hi2 = tmpvec2%tmpvec2;
        d = gold_rat3(n, ni, ni2, S, hi2, idmx,tol);
        oneplnid = (1+ ni * d[0]);
        temptcom.col(0) = sx.col(0)/oneplnid;
        temptcom.col(1) = sx.col(1)/oneplnid;
        tcom = -d[0] * temptcom.t();
        A = xx + tcom * sx;
        B = sxy + tcom * sy;
        down = A(0,0) * A(1,1) - A(0,1) * A(0,1);
        b2(0) = (A(1,1) * B(0) - A(0,1) * B(1))/down;
        b2(1) = (- A(0,1) * B(0) + A(0,0) * B(1))/down;
      }
      if(ret==1){
        se = (S - d[0] * sum(ni2 % hi2/ oneplnid ) )/n;
        seb = A(0,0) / down * se;
        stat = b2(1)*b2(1)/ seb;

        mymat(i,0) = stat;
        mymat(i,1) = R::pf(mymat(i,0), 1, n-4, false, logged);
      }
      else if(ret == 2){
        se = (S - d[0] * sum(ni2 % hi2/ oneplnid ) )/n;
        seb = A(0,0) / down * se;
        stat = b2(1)*b2(1)/ seb;

        info1 = -0.5 * d[1]-0.5*n*(1.837877-logn+1);
        info2 = -2 * info1;
        bic = info2 + (D + 2) * logn;

        mymat(i,0) = stat;
        mymat(i,1) = R::pf(mymat(i,0), 1, n, false, logged);
        mymat(i,2) = bic;
      }
      else{
        info1 = -0.5 * d[1]-0.5*n*(1.837877-logn+1);
        info2 = -2 * info1;
        bic = info2 + (D + 2) * logn;
        mymat(i,0) = bic;
      }
    }
    #ifdef _OPENMP
    }
    #endif
  }
  else{
    vec Xi(n), sxy(2), b1(2), tmpvec(n), tmpvec2(idmx), hi2(idmx),b2(2),B(2), mx(idmx);
    mat sx(idmx,2), temptcom(idmx,2), tcom(2,idmx), A(2,2);
    vec oneplnid;
    sx.col(0) = ni;
    sx.col(1) = ni;
    int ij=0;
    sxy[0] = Sy;

    double S=0,down=0, se=0,seb=0,info1=0,info2=0,stat,bic;
    vec d(2);
    mat  xx(2,2);
    xx(0,0) = n;

    for(int i = 0; i < D; i++) {
      Xi = x.col(i);
      xx(0,1) = xs[i];
      xx(1,0) = xs[i];
      xx(1,1) = xs2[i];
      sx.col(1) =  group_sum_helper<vec,vec,IntegerVector>(Xi, id, &idmn,&idmx);
      sxy[1] = sum(Xi % y);
      mx = sx.col(1)/ni;
      b1[0] = be.row(i)[0];
      b1[1] = be.row(i)[1];
      tmpvec = y - b1(0) - b1(1) * Xi;
      S = sum_with< square2<double>, vec>(tmpvec);
      tmpvec2 = my - b1(0) - b1(1) * mx;
      hi2 = tmpvec2%tmpvec2;
      d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
      oneplnid = (1+ ni * d[0]);
      temptcom.col(0) = sx.col(0)/oneplnid;
      temptcom.col(1) = sx.col(1)/oneplnid;
      tcom = -d[0] * temptcom.t();

      A = xx + tcom * sx;

      B = sxy + tcom * sy;

      down = A(0,0) * A(1,1) - A(0,1)*A(0,1);
      b2(0) = (A(1,1) * B(0) - A(0,1) * B(1))/down;
      b2(1) = (- A(0,1) * B(0) + A(0,0) * B(1))/down;
      ij = 2;
      while(ij++<maxiters && sum_with<std::abs,vec>(b2 - b1)>tol){
        b1 = b2;
        tmpvec = y - b1(0) - b1(1) * Xi;
        S = sum_with< square2<double>, vec>(tmpvec);
        tmpvec2 = my - b1(0) - b1(1) * mx;
        hi2 = tmpvec2%tmpvec2;
        d = gold_rat3(n, ni, ni2, S, hi2, idmx,tol);
        oneplnid = (1+ ni * d[0]);
        temptcom.col(0) = sx.col(0)/oneplnid;
        temptcom.col(1) = sx.col(1)/oneplnid;
        tcom = -d[0] * temptcom.t();
        A = xx + tcom * sx;
        B = sxy + tcom * sy;
        down = A(0,0) * A(1,1) - A(0,1) * A(0,1);
        b2(0) = (A(1,1) * B(0) - A(0,1) * B(1))/down;
        b2(1) = (- A(0,1) * B(0) + A(0,0) * B(1))/down;
      }
      if(ret==1){
        se = (S - d[0] * sum(ni2 % hi2/ oneplnid ) )/n;
        seb = A(0,0) / down * se;
        stat = b2(1)*b2(1)/ seb;

        mymat(i,0) = stat;
        mymat(i,1) = R::pf(mymat(i,0), 1, n-4, false, logged);
      }
      else if(ret == 2){
        se = (S - d[0] * sum(ni2 % hi2/ oneplnid ) )/n;
        seb = A(0,0) / down * se;
        stat = b2(1)*b2(1)/ seb;

        info1 = -0.5 * d[1]-0.5*n*(1.837877-logn+1);
        info2 = -2 * info1;
        bic = info2 + (D + 2) * logn;

        mymat(i,0) = stat;
        mymat(i,1) = R::pf(mymat(i,0), 1, n, false, logged);
        mymat(i,2) = bic;
      }
      else{
        info1 = -0.5 * d[1]-0.5*n*(1.837877-logn+1);
        info2 = -2 * info1;
        bic = info2 + (D + 2) * logn;
        mymat(i,0) = bic;
      }
    }
  }

  return mymat;
}

RcppExport SEXP Rfast_rint_regs(SEXP XSEXP,SEXP YSEXP,SEXP idSEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP parallelSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< NumericVector >::type Y(YSEXP);
  traits::input_parameter< IntegerVector >::type id(idSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const bool >::type logged(loggedSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = rint_regs(X,Y,id,tol,logged,parallel,maxiters);
  return __result;
  END_RCPP
}


////////////////////////////////////////////////////////////////////////////////////////


//[[Rcpp::export]]
List rint_mle(NumericVector X, IntegerVector id, const bool ranef, const double tol, const int maxiters){
  int n = X.size(),idmx,idmn;
  maximum<int>(id.begin(),id.end(),idmx);
  minimum<int>(id.begin(),id.end(),idmn);
  vec x(X.begin(),n,false);
  double sxy = sum(x);
  vec sx = group_sum_helper<vec,vec,IntegerVector>(x, id, &idmn,&idmx);

  List res;

  vec ni=Tabulate<vec,IntegerVector>(id,idmx);

  if (var(ni) == 0) {
    List tmp;

    tmp = varcomps_mle(X,id,idmx,tol);
    NumericVector mat = tmp["mat"];
    NumericVector info(3);
    info(0) = mat(0);
    info(1) = mat(1);
    info(2) = mat(2);
    double tmpd = mat(3);
    res["info"] = info;
    NumericVector syina = tmp["syina"];

    if(ranef)
      res["ranef"] = mat(0)/(mat(0) + mat(1)/tmpd) * syina/tmpd;
  }
  else {
    vec mx = sx/ni;
    vec com = ni % sx;
    double b1 = sxy/n;
    vec xminb1 = x - b1;
    double S = sum(xminb1%xminb1);
    vec mxminb1 = mx-b1;
    vec hi2 = mxminb1%mxminb1;

    vec d(2);
    vec ni2 = ni%ni;
    d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
    vec oneplnid = 1 + ni * d[0];
    double down = n - d[0] * sum(ni2/(oneplnid));
    double b2 = (sxy - d[0] * sum(com/(oneplnid)))/down;
    int i = 2;
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
    NumericVector info(3);

    double sigma = S/n;
    info[1] = sigma/(1 + d[0]);
    info[0] = sigma - info[1];
    info[2] = -0.5 * (d[1] + n * (1.83787706640935-log(n) + 1));
    res["info"] = info;
    res["my"] = b2;
    if (ranef) {
      res["ranef"] = conv_to<rowvec>::from((mx - b2) % (d[0] * ni/(oneplnid)));
    }
  }
  return res;
}

RcppExport SEXP Rfast_rint_mle(SEXP XSEXP,SEXP idSEXP,SEXP ranefSEXP,SEXP tolSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericVector >::type X(XSEXP);
  traits::input_parameter< IntegerVector >::type id(idSEXP);
  traits::input_parameter< const bool >::type ranef(ranefSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = rint_mle(X,id,ranef,tol,maxiters);
  return __result;
  END_RCPP
}
