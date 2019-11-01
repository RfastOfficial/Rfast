//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

RcppExport SEXP Rfast_transpose(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = wrap(Rfast::matrix::transpose(x));
    return __result;
END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////////


IntegerVector mat_mat(NumericMatrix x,NumericMatrix y){
  const int n=x.ncol(),p=y.ncol();
  LogicalMatrix f(p,n);
  NumericVector tmp;
  for(int i=0;i<n;++i){
    tmp=x.column(i);
    for(int j=0;j<p;++j){
      f(j,i)=as<bool>(all(tmp==y.column(j)));
    }
  }
  return colSums(f);
}

RcppExport SEXP Rfast_mat_mat(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericMatrix >::type y(ySEXP);
    __result = wrap(mat_mat(x,y));
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////////////////////////////////


RcppExport SEXP Rfast_mat_mult_p(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericMatrix >::type y(ySEXP);
    __result = wrap(Rfast::matrix::matrix_multiplication(x,y));
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////////////////////////////////


NumericMatrix submatrix(NumericMatrix x,const int rowstart,const int rowend,const int colstart,const int colend){
  return x(Range(rowstart-1,rowend-1),Range(colstart-1,colend-1));
}

RcppExport SEXP Rfast_submatrix(SEXP xSEXP,SEXP rowstartSEXP,SEXP rowendSEXP,SEXP colstartSEXP,SEXP colendSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const int >::type rowstart(rowstartSEXP);
    traits::input_parameter< const int >::type rowend(rowendSEXP);
    traits::input_parameter< const int >::type colstart(colstartSEXP);
    traits::input_parameter< const int >::type colend(colendSEXP);
    __result = wrap(submatrix(x,rowstart,rowend,colstart,colend));
    return __result;
END_RCPP
}


double sum_XopY(SEXP x,SEXP y,const char oper){
  switch(oper){
    case '+': return sum_x_op_y< madd<double>,madd<double> >(x,y);
    case '-': return sum_x_op_y< mdiff<double>,madd<double> >(x,y);
    case '*': return sum_x_op_y< mmult<double>,madd<double> >(x,y);
    case '/': return sum_x_op_y< mdiv<double>,madd<double> >(x,y);
    default: stop("The operation doesn't supported.");
  }
  return 0.0;
}

RcppExport SEXP Rfast_sum_XopY(SEXP x,SEXP y,SEXP operSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const char  >::type oper(operSEXP);
    __result = wrap(sum_XopY(x,y,oper));
    return __result;
END_RCPP
}

////////////////////////////////////////////////////////////////////////


double sum_XopX(SEXP x,const char oper){
  switch(oper){
    case '+': return sum_x_op_x< madd<double>,madd<double> >(x);
    case '-': return sum_x_op_x< mdiff<double>,madd<double> >(x);
    case '*': return sum_x_op_x< mmult<double>,madd<double> >(x);
    case '/': return sum_x_op_x< mdiv<double>,madd<double> >(x);
    default: stop("The operation doesn't supported.");
  }
  return 0.0;
}

RcppExport SEXP Rfast_sum_XopX(SEXP x,SEXP operSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const char  >::type oper(operSEXP);
    __result = wrap(sum_XopX(x,oper));
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////////////////