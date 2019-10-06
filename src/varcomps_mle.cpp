//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;

// [[Rcpp::export]]
List varcomps_mle(NumericVector x,IntegerVector ina,const int n,const double tol) {
  const double pi=3.14159265359;
  const int N=x.size(),d=N/n;
  NumericVector y=minus_mean(x,mean(x)),syina=group_sum_helper<NumericVector,NumericVector,IntegerVector>(y,ina,nullptr,nullptr);
  double sy2=sum_with< square2<double>,NumericVector>(syina),a=0,ratio=2.0/(sqrt(5) + 1),
  sy=sum_with< square2<double>,NumericVector>(y),b=sy/N,s=b;
  double x1=b-ratio*b,x2=ratio*b;
  double se=s-x1;
  double f1=N*log(se)+n*log1p(d*x1/se)+sy/se-x1/(se*se+d*x1*se)*sy2; 
  se=s-x2;
  double f2=N*log(se)+n*log1p(d*x2/se)+sy/se-x2/(se*se+d*x2*se)*sy2;
  
  while (abs(b-a)>tol){
    if(f2>f1){
      b=x2;
      x2=x1;
      f2=f1;
      x1=b - ratio * (b - a);
      se=s - x1;
      f1=N * log(se) + n * log1p(d * x1 / se) + sy/se - x1 / (se*se + d * x1 * se) * sy2 ;
    } else {
      a=x1;
      x1=x2;
      f1=f2;
      x2=a + ratio * (b - a);
      se=s - x2;
      f2=N * log(se) + n * log1p(d * x2 / se) + sy/se - x2 / (se*se + d * x2 * se) * sy2; 
    }
  }
  const double tau=(a+b)/2.0;
  NumericVector m(4);
  m[0]=tau;
  m[1]=s - tau;
  m[2]=-0.5 * f2 - N*0.5 * log(2 * pi);
  m[3]=d;
  List f;
  f["syina"]=syina;
  f["mat"]=m;
  return f;
}

RcppExport SEXP Rfast_varcomps_mle(SEXP xSEXP,SEXP inaSEXP,SEXP nSEXP,SEXP tolSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< IntegerVector >::type ina(inaSEXP);
    traits::input_parameter< const int >::type n(nSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    __result = wrap(varcomps_mle(x,ina,n,tol));
    return __result;
END_RCPP
}