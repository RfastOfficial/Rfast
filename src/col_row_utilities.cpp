
//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"
#include <chrono>
#include <random>
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
SEXP col_all(SEXP x){
  const int n=Rf_ncols(x),p=Rf_nrows(x);
  SEXP f=Rf_allocVector(LGLSXP,n);
  int *start=LOGICAL(x),*end=start+p,*ff=LOGICAL(f);
  for(int i=0;i<n;++i,++ff){
    *ff=my_all(start,end);
    start=end;
    end+=p;
  }
  return f;
}

RcppExport SEXP Rfast_col_all(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_all(x);
    return __result;
END_RCPP
}

LogicalVector row_all(LogicalMatrix x){
  const int n=x.nrow();
  LogicalVector f(n);
  for(int i=0;i<n;++i)
    f[i]=as<bool>(all(x.row(i)));
  return f;
}

RcppExport SEXP Rfast_row_all(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< LogicalMatrix >::type x(xSEXP);
    __result = row_all(x);
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////////



SEXP col_any(SEXP x){
  const int n=Rf_ncols(x),p=Rf_nrows(x);
  SEXP f=Rf_allocVector(LGLSXP,n);
  int *start=LOGICAL(x),*end=start+p,*ff=LOGICAL(f);
  for(int i=0;i<n;++i,++ff){
    *ff=my_any(start,end);
    start=end;
    end+=p;
  }
  return f;
}

SEXP row_any(SEXP x){
  int nrow=Rf_nrows(x);
  SEXP F=PROTECT(Rf_allocVector(LGLSXP,nrow));
  int *xx=INTEGER(x),*endx=xx+LENGTH(x),*f=INTEGER(F),*startx,*startf;
  const int *endf=f+LENGTH(F);
  for(startf=f;startf!=endf;++startf)
    *startf=0;
  while(xx!=endx){
    for(startf=f,startx=xx,xx+=nrow;startx!=xx;++startf,++startx){
      if(*startx){
        *startf=1;
      }
    }
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_row_any(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = row_any(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_col_any(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_any(x);
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////


IntegerVector col_count_values(NumericMatrix x,NumericVector values){
  const int n=values.size();
  IntegerVector f(n);
  for(int i=0;i<n;++i){
    f[i]=count_value_helper<NumericVector,double>(x.column(i),values[i]);
  }
  return f;
}

RcppExport SEXP Rfast_col_count_values(SEXP xSEXP,SEXP valuesSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type values(valuesSEXP);
    __result = col_count_values(x,values);
    return __result;
END_RCPP
}


IntegerVector row_count_values(NumericMatrix x,NumericVector values){
  const int n=values.size();
  IntegerVector f(n);
  for(int i=0;i<n;++i){
    f[i]=count_value_helper<NumericVector,double>(x.row(i),values[i]);
  }
  return f;
}

RcppExport SEXP Rfast_row_count_values(SEXP xSEXP,SEXP valuesSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type values(valuesSEXP);
    __result = row_count_values(x,values);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////


SEXP col_false(SEXP x){
  const int p=Rf_ncols(x),n=Rf_nrows(x);
  SEXP f=Rf_allocVector(INTSXP,p);
  int *ff=INTEGER(f),*xx=LOGICAL(x),*endx=xx+LENGTH(x);
  for(;xx!=endx;xx+=p,++ff)
    *ff=n-True(xx,xx+p);
  return f;
}

SEXP row_false(SEXP x){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F=PROTECT(Rf_allocVector(INTSXP,nrow));
  int *xx=INTEGER(x),*end=xx+ncol*nrow,*f=INTEGER(F),*startx,*startf;
  const int *endf=f+LENGTH(F);
  for(startf=f;startf!=endf;++startf)
    *startf=ncol;
  while(xx!=end){
    for(startf=f,startx=xx,xx+=nrow;startx!=xx;++startf,++startx){
      *startf-=*startx;
    }
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_row_false(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = row_false(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_col_false(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_false(x);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////


IntegerVector col_len_sort_un_int(IntegerMatrix x){
  const int p=x.ncol();
  IntegerVector f(p);
  for(int i=0;i<p;++i)
    f[i]=len_sort_unique_int(x.column(i));
  return f;
}

IntegerVector row_len_sort_un_int(IntegerMatrix x){
  const unsigned int p=x.nrow();
  IntegerVector F(p);
  IntegerVector::iterator FF=F.begin();
    for(int i=0;FF!=F.end();++FF,++i){
    *FF=len_sort_unique_int(x.row(i));
  }
  return F;
}

RcppExport SEXP Rfast_row_len_sort_un_int(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    __result = row_len_sort_un_int(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_col_len_sort_un_int(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    __result = col_len_sort_un_int(x);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////


RcppExport SEXP Rfast_col_mads(SEXP xSEXP,SEXP methodSEXP,SEXP na_rmSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const string >::type method(methodSEXP);
    traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = Rfast::matrix::colMads(x,method,na_rm,parallel);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_row_mads(SEXP xSEXP,SEXP methodSEXP,SEXP na_rmSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const string >::type method(methodSEXP);
    traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = Rfast::matrix::rowMads(x,method,na_rm,parallel);
    return __result;
END_RCPP
}


////////////////////////////////////////////////////


SEXP col_max_indices(NumericMatrix x){
  unsigned int i=0,p=x.ncol();
  arma::mat X = arma::mat(x.begin(), x.nrow(), p, false);
  SEXP F=PROTECT(Rf_allocVector(INTSXP,p));
  int *FF=INTEGER(F);
  for(;i<p;++i,++FF)
    *FF=(X.col(i)).index_max()+1;
  UNPROTECT(1);
  return F;
}

SEXP col_max(SEXP x){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F;
  switch(TYPEOF(x)){
    case REALSXP:{
      F=PROTECT(Rf_allocVector(REALSXP,ncol));
      double *xx=REAL(x),*end=xx+ncol*nrow,*f=REAL(F);
      for(;xx!=end;xx+=nrow,++f)
        maximum<double>(xx,xx+nrow,*f);
      break;
    }
    default:{
      F=PROTECT(Rf_allocVector(INTSXP,ncol));
      int *xx=INTEGER(x),*end=xx+ncol*nrow,*f=INTEGER(F);
      for(;xx!=end;xx+=nrow,++f)
        maximum<int>(xx,xx+nrow,*f);
      break;
    }
  }
  UNPROTECT(1);
  return F;
}

// find the maximum index of its collumn
RcppExport SEXP Rfast_col_max_indices(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = col_max_indices(x);
    return __result;
END_RCPP
}

// find the maximum value of its collumn
RcppExport SEXP Rfast_col_max(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_max(x);
    return __result;
END_RCPP
}



// [[Rcpp::export]]
SEXP row_max_indices(NumericMatrix x){
  const int p=x.nrow();
  mat X = mat(x.begin(), p, x.ncol(), false); 
  SEXP F=PROTECT(Rf_allocVector(INTSXP,p));
  int *FF=INTEGER(F);
  for(int i=0;i<p;++i,++FF)
    *FF=(X.row(i)).index_max()+1;
  UNPROTECT(1);
  return F;
}

SEXP row_max(SEXP x){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F;
  F=PROTECT(Rf_allocVector(REALSXP,nrow));
  double *xx=REAL(x),*end=xx+ncol*nrow,*f=REAL(F),*x3,*ff;
  const double *endf=f+LENGTH(F);
  for(ff=f;ff!=endf;++ff,++xx)
    *ff=*xx;
  for(;xx!=end;)
    for(ff=f,x3=xx,xx+=nrow;x3!=xx;++ff,++x3){
      *ff=std::max(*ff,*x3);
    }
    UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_row_max_indices(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = row_max_indices(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_row_max(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = row_max(x);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////////


arma::rowvec col_means(NumericMatrix x){
  mat X = mat(x.begin(), x.nrow(), x.ncol(), false); 
  return mean(X, 0); 
}

colvec row_means(NumericMatrix x){
  mat X = mat(x.begin(), x.nrow(), x.ncol(), false); 
  return mean(X, 1); 
}

RcppExport SEXP Rfast_row_means(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = row_means(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_col_means(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = col_means(x);
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////


// colMedians
RcppExport SEXP Rfast_col_meds(SEXP xSEXP,SEXP na_rmSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = Rfast::matrix::colMedian(x,na_rm,parallel);
    return __result;
END_RCPP
}

// rowMedians
RcppExport SEXP Rfast_row_meds(SEXP xSEXP,SEXP na_rmSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = Rfast::matrix::rowMedian(x,na_rm,parallel);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////


SEXP col_min_indices(NumericMatrix x){
  unsigned int i=0,p=x.ncol();
  mat X = mat(x.begin(), x.nrow(), p, false);
  SEXP F=PROTECT(Rf_allocVector(INTSXP,p));
  int *FF=INTEGER(F);
  for(;i<p;++i,++FF)
    *FF=(X.col(i)).index_min()+1;
  UNPROTECT(1);
  return F;
}

SEXP col_min(SEXP x){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F;
  switch(TYPEOF(x)){
    case REALSXP:{
      F=PROTECT(Rf_allocVector(REALSXP,ncol));
      double *xx=REAL(x),*end=xx+ncol*nrow,*f=REAL(F);
    for(;xx!=end;xx+=nrow,++f)
        minimum<double>(xx,xx+nrow,*f);
      break;
    }
    default:{
      F=PROTECT(Rf_allocVector(INTSXP,ncol));
    int *xx=INTEGER(x),*end=xx+ncol*nrow,*f=INTEGER(F);
    for(;xx!=end;xx+=nrow,++f)
        minimum<int>(xx,xx+nrow,*f);
      break;
    }
  }
  UNPROTECT(1);
  return F;
}

// find the minimum index of its collumn
RcppExport SEXP Rfast_col_min_indices(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = col_min_indices(x);
    return __result;
END_RCPP
}

// find the minimum value of its collumn
RcppExport SEXP Rfast_col_min(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_min(x);
    return __result;
END_RCPP
}



// [[Rcpp::export]]
NumericVector row_min_indices(NumericMatrix x){
  const unsigned int p=x.nrow();
  mat X = mat(x.begin(), p, x.ncol(), false); 
  NumericVector F(p);
  NumericVector::iterator FF=F.begin();
  for(unsigned int i=0;i<p;++i,++FF)
      *FF=(X.row(i)).index_min()+1;
  return F;
}

SEXP row_min(SEXP x){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F;
  F=PROTECT(Rf_allocVector(REALSXP,nrow));
  double *xx=REAL(x),*end=xx+ncol*nrow,*f=REAL(F),*x3,*ff;
  const double *endf=f+LENGTH(F);
  for(ff=f;ff!=endf;++ff,++xx)
    *ff=*xx;
  for(;xx!=end;)
    for(ff=f,x3=xx,xx+=nrow;x3!=xx;++ff,++x3){
      *ff=std::min(*ff,*x3);
    }
    UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_row_min_indices(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = row_min_indices(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_row_min(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = row_min(x);
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////


SEXP col_min_max(SEXP x){
  const int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F;
  switch(TYPEOF(x)){
    case REALSXP:{
      F=PROTECT(Rf_allocMatrix(REALSXP,2,ncol));
      double *xx=REAL(x),*end=xx+LENGTH(x),*f=REAL(F),min,max;
      for(;xx!=end;xx+=nrow,f+=2){
        min_max<double>(xx,xx+nrow,min,max);
      *f=min;
      f[1]=max;
      }
      break;
    }
    default:{
      F=PROTECT(Rf_allocMatrix(INTSXP,2,ncol));
      int *xx=INTEGER(x),*end=xx+LENGTH(x),*f=INTEGER(F),min,max;
    for(;xx!=end;xx+=nrow,f+=2){
        min_max<int>(xx,xx+nrow,min,max);
        *f=min;
        f[1]=max;
    }
    }
  }
  UNPROTECT(1);
  return F;
}


RcppExport SEXP Rfast_col_min_max(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_min_max(x);
    return __result;
END_RCPP
}



SEXP row_min_max(SEXP x){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F;
  F=PROTECT(Rf_allocMatrix(REALSXP,2,nrow));
  double *xx=REAL(x),*end=xx+ncol*nrow,*f=REAL(F),*x3,*ff;
  const double *endf=f+(nrow<<1);
  for(ff=f;ff!=endf;ff+=2,++xx)
    *ff=ff[1]=*xx;
  for(;xx!=end;)
    for(ff=f,x3=xx,xx+=nrow;x3!=xx;ff+=2,++x3){
      if(*ff>*x3)
        *ff=*x3;
      else if(ff[1]<*x3)
        ff[1]=*x3;
    }
    UNPROTECT(1);
  return F;
}


RcppExport SEXP Rfast_row_min_max(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = row_min_max(x);
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////////


using std::nth_element;

SEXP col_nth(NumericMatrix x,IntegerVector elems,const bool descend,const bool na_rm,const bool index){
    const int n=elems.size();
    SEXP F;
    NumericVector y(x.nrow());
    IntegerVector::iterator els=elems.begin();
    if(index){
      F=PROTECT(Rf_allocVector(INTSXP,n));
      int *ff=INTEGER(F);
      for(int i=0;i!=n;++ff,++i,++els){
          y=x.column(i);
          *ff=nth_helper_index<NumericVector>(y,*els,descend,na_rm);
      }
    }else{
      F=PROTECT(Rf_allocVector(REALSXP,n));
      double *ff=REAL(F);
      for(int i=0;i!=n;++ff,++i,++els){
          y=x.column(i);
          *ff=nth_helper<NumericVector>(y,*els,descend,na_rm);
      }
    }
    UNPROTECT(1);
    return F;
}

// nth_element
RcppExport SEXP Rfast_col_nth(SEXP xSEXP,SEXP ySEXP,SEXP descendSEXP,SEXP na_rmSEXP,SEXP indexSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< IntegerVector >::type y(ySEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    traits::input_parameter< const bool >::type index(indexSEXP);
    __result = col_nth(x,y,descend,na_rm,index);
    return __result;
END_RCPP
}

SEXP row_nth(NumericMatrix x,IntegerVector elems,const bool descend,const bool na_rm,const bool index){
    const int n=elems.size();
    NumericVector y(x.ncol());
    SEXP F;
    IntegerVector::iterator els=elems.begin();
    if(index){
      F=PROTECT(Rf_allocVector(INTSXP,n));
      int *ff=INTEGER(F);
      for(int i=0;i!=n;++ff,++i,++els){
        y=x.row(i);
        *ff=nth_helper_index<NumericVector>(y,*els,descend,na_rm);
      }
    }else{
      F=PROTECT(Rf_allocVector(REALSXP,n));
      double *ff=REAL(F);
      for(int i=0;i!=n;++ff,++i,++els){
        y=x.row(i);
        *ff=nth_helper<NumericVector>(y,*els,descend,na_rm);
      }
    }
    UNPROTECT(1);
    return F;
}

RcppExport SEXP Rfast_row_nth(SEXP xSEXP,SEXP ySEXP,SEXP descendSEXP,SEXP na_rmSEXP,SEXP indexSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< IntegerVector >::type y(ySEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    traits::input_parameter< const bool >::type index(indexSEXP);
    __result = row_nth(x,y,descend,na_rm,index);
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////////////////


IntegerMatrix col_order(NumericMatrix x,const bool stable,const bool descending){
    const int ncl=x.ncol();
    IntegerMatrix f(x.nrow(),ncl);
    for(int i=0;i<ncl;++i){
      f.column(i)=Order(x.column(i),stable,descending);
    }
    return f;
}

RcppExport SEXP Rfast_col_order(SEXP xSEXP,SEXP stableSEXP,SEXP descendingSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type stable(stableSEXP);
    traits::input_parameter< const bool >::type descending(descendingSEXP);
    __result = col_order(x,stable,descending);
    return __result;
END_RCPP
}

IntegerMatrix row_order(NumericMatrix x,const bool stable,const bool descending){
  const int nrw=x.nrow();
  IntegerMatrix f(nrw,x.ncol());
    for(int i=0;i<nrw;++i)
      f.row(i)=Order(x.row(i),stable,descending);
  return f;
}

RcppExport SEXP Rfast_row_order(SEXP xSEXP,SEXP stableSEXP,SEXP descendingSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type stable(stableSEXP);
    traits::input_parameter< const bool >::type descending(descendingSEXP);
    __result = row_order(x,stable,descending);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////////////////////////////


NumericVector col_prods(SEXP x,string method){
	const int n=Rf_ncols(x);
	NumericVector f(n);
	if(method == "direct"){
	  mat X(REAL(x), Rf_nrows(x), n, false);
	  rowvec ff(f.begin(),n,false);
	  ff = prod(X, 0);
	}else if(method =="expsumlog"){
		int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
    	double *xx=REAL(x),*end=xx+ncol*nrow,*ff=f.begin();
	    for(;xx!=end;xx+=nrow,++ff){
	        *ff=exp(accumulate(xx,xx+nrow,0.0,[](double& s,double x){return x<0 ? s+x : s+log(x);}));
	    }
	}else{
		stop("Error: Unsupported method.");
	}
  	return f;
}

RcppExport SEXP Rfast_col_prods(SEXP x,SEXP methodSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< string >::type method(methodSEXP);
    __result = col_prods(x,method);
    return __result;
END_RCPP
}


NumericVector row_prods(NumericMatrix x){
    const int n=x.nrow();
	NumericVector f(n);
  	mat X = mat(x.begin(), n, x.ncol(), false); 
	colvec ff(f.begin(),n,false);
	ff = prod(X, 1);
  	return f; 
}

RcppExport SEXP Rfast_row_prods(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = row_prods(x);
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////////////


using std::string;


NumericMatrix col_ranks(NumericMatrix x,string method,const bool descend,const bool stable){
  const int n=x.ncol();
  NumericMatrix f(x.nrow(),n);
  for(int i=0;i<n;++i){
    f.column(i)=Rank(x.column(i),method,descend,stable);
  }
  return f; 
}

RcppExport SEXP Rfast_col_ranks(SEXP xSEXP,SEXP methodSEXP,SEXP descendSEXP,SEXP stableSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< string >::type method(methodSEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    traits::input_parameter< const bool >::type stable(stableSEXP);    
    __result = col_ranks(x,method,descend,stable);
    return __result;
END_RCPP
}

NumericMatrix row_ranks(NumericMatrix x,string method,const bool descend,const bool stable){
  const int n=x.nrow();
  NumericMatrix f(n,x.ncol());
  for(int i=0;i<n;++i){
    f.row(i)=Rank(x.row(i),method,descend,stable);
  }
  return f; 
}

RcppExport SEXP Rfast_row_ranks(SEXP xSEXP,SEXP methodSEXP,SEXP descendSEXP,SEXP stableSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< string >::type method(methodSEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    traits::input_parameter< const bool >::type stable(stableSEXP);
    __result = row_ranks(x,method,descend,stable);
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////////////////////

RcppExport SEXP Rfast_row_shuffle(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rfast::matrix::rowShuffle(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_col_shuffle(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rfast::matrix::colShuffle(x);
    return __result;
END_RCPP
}

////////////////////////////////////////////////////////////////


template<class Ret,class T1,class F1,class F2>
Ret col_sums(T1 x,SEXP indices){
    const int n=Rf_isNull(indices) ? 0 : LENGTH(indices);
    F1 X(x.begin(), x.nrow(), x.ncol(), false);
    Ret f(n==0 ? X.n_cols : n);
    if(n==0){
        F2 ff(f.begin(),X.n_cols,false);
        ff = sum(X,0);
    }else{
        IntegerVector ind(indices);
        for(int i=0;i<n;++i)
            f[i]=accu(X.col(ind[i]-1));
    }
    return f;
}

template<class Ret,class T1,class F1,class F2>
Ret row_sums(T1 x,SEXP indices){
  const int n=Rf_isNull(indices) ? 0 : LENGTH(indices);
  F1 X(x.begin(), x.nrow(), x.ncol(), false);
  Ret f(n==0 ? X.n_rows : n);
  if(n==0){
    F2 ff(f.begin(),X.n_rows,false);
    ff = sum(X,1);
  }else{
    IntegerVector ind(indices);
    for(int i=0;i<n;++i)
      f[i]=accu(X.row(ind[i]-1));
  }
  return f;
}

RcppExport SEXP Rfast_row_sums(SEXP x,SEXP indices) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = Rf_isInteger(x) ? row_sums<IntegerVector,IntegerMatrix,imat,icolvec>(x,indices)
        : row_sums<NumericVector,NumericMatrix,mat,colvec>(x,indices);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_col_sums(SEXP x,SEXP indices) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = Rf_isInteger(x) ? col_sums<IntegerVector,IntegerMatrix,imat,irowvec>(x,indices)
    		   : col_sums<NumericVector,NumericMatrix,mat,rowvec>(x,indices);
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////


IntegerMatrix col_tabulate(IntegerMatrix x,int nroww){
  const int ncl=x.ncol();
  IntegerMatrix f(nroww,ncl);
  for(int i=0;i<ncl;++i)
    f.column(i)=Tabulate<IntegerVector>(x.column(i),nroww);
  return f;
}

IntegerMatrix row_tabulate(IntegerMatrix x,int ncoll){
  const int nrw=x.nrow();
  IntegerMatrix f(nrw,ncoll);
  for(int i=0;i<nrw;++i)
    f.row(i)=Tabulate<IntegerVector>(x.row(i),ncoll);
  return f;
}

RcppExport SEXP Rfast_row_tabulate(SEXP xSEXP,SEXP ncollSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    traits::input_parameter< int >::type ncoll(ncollSEXP);
    __result = row_tabulate(x,ncoll);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_col_tabulate(SEXP xSEXP,SEXP nrowwSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    traits::input_parameter< int >::type nroww(nrowwSEXP);
    __result = col_tabulate(x,nroww);
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////


SEXP col_true(SEXP x){
  const int p=Rf_nrows(x);
  SEXP f=Rf_allocVector(INTSXP,p);
  int *ff=INTEGER(f),*xx=LOGICAL(x),*endx=xx+LENGTH(x);
  for(;xx!=endx;xx+=p,++ff)
    *ff=True(xx,xx+p);
  return f;
}

SEXP row_true(SEXP x){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F=PROTECT(Rf_allocVector(INTSXP,nrow));
  int *xx=INTEGER(x),*end=xx+ncol*nrow,*f=INTEGER(F),*startx,*startf;
  const int *endf=f+LENGTH(F);
  for(startf=f;startf!=endf;++startf)
    *startf=0;
  while(xx!=end){
    for(startf=f,startx=xx,xx+=nrow;startx!=xx;++startf,++startx){
      *startf+=*startx;
    }
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_row_true(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = row_true(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_col_true(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_true(x);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////


SEXP col_true_false(SEXP x){
  const int n=Rf_nrows(x);
  SEXP f=Rf_allocMatrix(INTSXP,2,Rf_ncols(x));
  int *ff=INTEGER(f),*xx=LOGICAL(x),*endx=xx+LENGTH(x),t;
  for(;xx!=endx;xx+=n,ff+=2){
    t=True(xx,xx+n);
    *ff=n-t;
    ff[1]=t;
  }
    return f;
}

RcppExport SEXP Rfast_col_true_false(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_true_false(x);
    return __result;
END_RCPP
}


SEXP row_true_false(SEXP x){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP F=PROTECT(Rf_allocMatrix(INTSXP,2,nrow));
  int *xx=INTEGER(x),*end=xx+ncol*nrow,*f=INTEGER(F),*startx,*startf;
  const int *endf=f+LENGTH(F);
  for(startf=f;startf!=endf;startf+=2){
    *startf=ncol;
    startf[1]=0;
  }
  while(xx!=end){
    for(startf=f,startx=xx,xx+=nrow;startx!=xx;startf+=2,++startx){
      *startf-=*startx;
      startf[1]+=*startx;
    }
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_row_true_false(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = row_true_false(x);
    return __result;
END_RCPP
}


////////////////////////////////////////////////

SEXP col_pmin(SEXP x,SEXP y){
  const int nrows=Rf_nrows(x),ncols=Rf_ncols(x);
  SEXP f=Rf_allocMatrix(REALSXP,nrows,ncols);
  double *startx=REAL(x),*end=startx+ncols*nrows,*starty=REAL(y),*startf=REAL(f),*endx;
  for(;startx!=end;){
    endx=startx+nrows;
    for(;startx!=endx;++startx,++starty,++startf)
      *startf=std::min(*startx,*starty);
  }
  return f;
}

RcppExport SEXP Rfast_col_pmin(SEXP x,SEXP y) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_pmin(x,y);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////

SEXP col_pmax(SEXP x,SEXP y){
  const int nrows=Rf_nrows(x),ncols=Rf_ncols(x);
  SEXP f=Rf_allocMatrix(REALSXP,nrows,ncols);
  double *startx=REAL(x),*end=startx+ncols*nrows,*starty=REAL(y),*startf=REAL(f),*endx;
  for(;startx!=end;){
    endx=startx+nrows;
    for(;startx!=endx;++startx,++starty,++startf)
      *startf=std::max(*startx,*starty);
  }
  return f;
}

RcppExport SEXP Rfast_col_pmax(SEXP x,SEXP y) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_pmax(x,y);
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////////

SEXP col_cum_prods(SEXP x){
  const int p=Rf_nrows(x);
  SEXP f=Rf_duplicate(x);
  double *ff=REAL(f),*endf=ff+LENGTH(f);
  int i=1;
  for(++ff;ff!=endf;++ff,++i){
    i!=p ? *ff*=ff[-1] : i=0;
  }
  return f;
}

RcppExport SEXP Rfast_col_cum_prods(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_cum_prods(x);
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////////

SEXP col_cum_maxs(SEXP x){
  int p=Rf_nrows(x);
  SEXP f=Rf_duplicate(x);
  double *ff=REAL(f),*endf=ff+LENGTH(f);
  int i=1;
  for(++ff;ff!=endf;++ff,++i){
    i!=p ? *ff=std::max(*ff,ff[-1]) : i=0;
  }
  return f;
}

RcppExport SEXP Rfast_col_cum_maxs(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_cum_maxs(x);
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////////

SEXP col_cum_sums(SEXP x){
  int p=Rf_nrows(x);
  SEXP f=Rf_duplicate(x);
  double *ff=REAL(f),*endf=ff+LENGTH(f);
  int i=1;
  for(++ff;ff!=endf;++ff,++i){
    i!=p ? *ff+=ff[-1] : i=0;
  }
  return f;
}

RcppExport SEXP Rfast_col_cum_sums(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_cum_sums(x);
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////////

SEXP col_cum_mins(SEXP x){
  int p=Rf_nrows(x);
  SEXP f=Rf_duplicate(x);
  double *ff=REAL(f),*endf=ff+LENGTH(f);
  int i=1;
  for(++ff;ff!=endf;++ff,++i){
    i!=p ? *ff=std::min(*ff,ff[-1]) : i=0;
  }
  return f;
}

RcppExport SEXP Rfast_col_cum_mins(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_cum_mins(x);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////////////////////////////


bool col_row_value(NumericMatrix X,double v){
  int i,n=X.nrow(),p=X.ncol();
  mat x(X.begin(),n,p,false);
  for(i=0;i<p;++i){
    if(all(x.col(i)==v)){
      return true;
    }
  }
  for(i=0;i<n;++i){
    if(all(x.row(i)==v)){
      return true;
    }
  }
  return false;
}

RcppExport SEXP Rfast_col_row_value(SEXP xSEXP,SEXP vSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const double >::type v(vSEXP);
    __result = col_row_value(x,v);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////////////////////////


NumericMatrix columns(NumericMatrix x,IntegerVector ind){
  const int nrw=x.nrow(),ncl=ind.size();
  NumericMatrix f(nrw,ncl);
  mat ff(f.begin(),nrw,ncl,false),xx(x.begin(),nrw,x.ncol(),false);
  for(int i=0;i<ncl;++i)
    ff.col(i)=xx.col(ind[i]-1);
  return f;
}

SEXP rows(SEXP X,SEXP Ind){
  const int nrw=Rf_nrows(X),ncl=Rf_ncols(X);
  SEXP F=PROTECT(Rf_allocMatrix(REALSXP,LENGTH(Ind),ncl));
  double *start = REAL(X),*ff=REAL(F),*xx=start;
  int *start_ind=INTEGER(Ind),*ind,*end_ind=start_ind+LENGTH(Ind);
  for(int i=0;i<ncl;++i){
    for(ind=start_ind;ind!=end_ind;++ind){
      xx=start+ *ind-1;
      *ff++=*xx;
    }
    start+=nrw;
  }
  UNPROTECT(1);
  return F;
}

RcppExport SEXP Rfast_columns(SEXP xSEXP,SEXP indSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< IntegerVector >::type ind(indSEXP);
    __result = columns(x,ind);
    return __result;
END_RCPP
}


RcppExport SEXP Rfast_rows(SEXP x,SEXP ind) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = rows(x,ind);
    return __result;
END_RCPP
}


RcppExport SEXP Rfast_col_vars(SEXP xSEXP,SEXP stdSEXP,SEXP na_rmSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type std(stdSEXP);
    traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = Rfast::matrix::colVars(x,std,na_rm,parallel);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_row_vars(SEXP xSEXP,SEXP stdSEXP,SEXP na_rmSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type std(stdSEXP);
    traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = Rfast::matrix::rowVars(x,std,na_rm,parallel);
    return __result;
END_RCPP
}
