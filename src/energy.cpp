
//Author: Manos Papadakis
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include "mn.h"

using namespace Rcpp;
using namespace arma;

List dcor(NumericMatrix x,NumericMatrix y) {
  NumericMatrix a = euclidean_dist(x,false);
  NumericMatrix b = euclidean_dist(y,false);
  const int n=a.ncol();
  mat aa(a.begin(),n,n,false);
  mat bb(b.begin(),n,n,false);
  rowvec ma = mean(aa,0);
  rowvec mb = mean(bb,0);
  mat A = aa.each_row() - ma;   
  A = A.each_col() - ma.t();   
  A = A + mean(ma);
  mat B = bb.each_row() - mb;  
  B = B.each_col() - mb.t();
  B = B + mean(mb);
  const double dcov = sqrt( mean(mean(A % B)) );
  const double dvarX = sqrt( mean(mean(square(A))) );
  const double dvarY = sqrt( mean(mean(square(B))) );
  const double dcor = dcov/sqrt(dvarX * dvarY);
  List l;
  l["dcov"] = dcov;
  l["dvarX"] = dvarX;
  l["dvarY"] = dvarY;
  l["dcor"] = dcor;
  return l;
}

RcppExport SEXP Rfast_dcor(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericMatrix >::type y(ySEXP);
    __result = wrap(dcor(x,y));
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////////////////////////


double dcov(NumericMatrix x,NumericMatrix y) {
  NumericMatrix a = euclidean_dist(x,false);
  NumericMatrix b = euclidean_dist(y,false);
  const int ncla=a.ncol(),nrwa=a.nrow();
  mat aa(a.begin(),nrwa,ncla,false);
  const int nclb=b.ncol(),nrwb=b.nrow();
  mat bb(b.begin(),nrwb,nclb,false);
  rowvec ma = mean(aa,0);
  rowvec mb = mean(bb,0);
  mat A = aa.each_row() - ma;   
  A = A.each_col() - ma.t();   
  A = A + mean(ma);
  mat B = bb.each_row() - mb;  
  B = B.each_col() - mb.t();
  B = B + mean(mb);
  return sqrt(mean(mean(A % B)));
}

RcppExport SEXP Rfast_dcov(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericMatrix >::type y(ySEXP);
    __result = dcov(x,y);
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////////////


double dvar(NumericMatrix x) {
  NumericMatrix a = euclidean_dist(x,false);
  const int ncla=a.ncol(),nrwa=a.nrow();
  mat aa(a.begin(),nrwa,ncla,false);
  rowvec ma = mean(aa,0);
  mat A = aa.each_row() - ma;   
  A = A.each_col() - ma.t();   
  A = A + mean(ma);
  return sqrt(mean(mean(square(A))));
}

RcppExport SEXP Rfast_dvar(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = dvar(x);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////////


double edist(NumericMatrix x,NumericMatrix y){
	const int n1=x.ncol(),n2=y.ncol();
	double mij=total_dista(x, y,false),mii=total_euclidean_dist(x,false),mjj=total_euclidean_dist(y,false);
	return (2 * mij - n2 * mii / n1 - n1 * mjj/n2 ) / (n1 + n2);
}

RcppExport SEXP Rfast_edist(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericMatrix >::type y(ySEXP);
    __result = edist(x,y);
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////////


static vec lower_tri(mat &x){
  int n=x.n_cols,i,j;
  vec f(n*(n-1)*0.5);
  vec::iterator ff=f.begin();
  for(i=0;i<n;++i)
    for(j=i+1;j<n;++j,++ff)
      *ff=x(j,i);
  return f;
}


double bcdcor(NumericMatrix x,NumericMatrix y) {
  NumericMatrix a = euclidean_dist(x,false);
  NumericMatrix b = euclidean_dist(y,false);
  const double n=a.ncol();
  mat aa(a.begin(),n,n,false);
  mat bb(b.begin(),n,n,false);
  rowvec ma = mean(aa,0),dgA(n);
  rowvec mb = mean(bb,0),dgB(n);
  const double mean_ma=mean(ma),mean_mb=mean(mb);
  double n_1=n/(n-1),n_2=n/(n-2);
  dgA=ma-mean_ma;
  dgB=mb-mean_mb;
  mat A = aa.each_row() - ma;   
  A = A.each_col() - ma.t();   
  A = A + mean_ma - aa/n;
  vec Al=lower_tri(A);
  mat B = bb.each_row() - mb;  
  B = B.each_col() - mb.t();
  B = B + mean_mb - bb/n;
  vec Bl=lower_tri(B);
  n_1*=n_1;
  n_2*=n_1;
  const double sdgab=sum(dgA%dgB),sdga=sum(square(dgA)),sdgb=sum(square(dgB));
  const double XY = n_1*(sum(Al%Bl)*2+sdgab) - n_2 * sdgab;
  const double XX = n_1*(sum(square(Al))*2+sdga) - n_2 * sdga;
  const double YY = n_1*(sum(square(Bl))*2+sdgb) - n_2 * sdgb;
  return XY / sqrt(XX * YY);
}

RcppExport SEXP Rfast_bcdcor(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericMatrix >::type y(ySEXP);
    __result = bcdcor(x,y);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////////////////////////////////