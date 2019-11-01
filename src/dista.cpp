//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "Rfast.h"
#include "mn.h"

using namespace arma;
using namespace Rcpp;
using std::string;

//[[Rcpp::export]]
NumericMatrix dista(NumericMatrix Xnew, NumericMatrix X,const bool sqr,const string type) {
  const int n=X.ncol(),nu=Xnew.ncol();
  mat xnew(Xnew.begin(),Xnew.nrow(),nu,false),x(X.begin(),X.nrow(),n,false);
  NumericMatrix disaa(n,nu);
  mat disa(disaa.begin(),n,nu,false);
  if(type == "euclidean"){
    if(sqr){
      for(int i=0;i<nu;++i){
        disa.col(i)=sum(square(x.each_col() - xnew.col(i)),0).t();
      }
    }else{
      for(int i=0;i<nu;++i){
        disa.col(i)=foreach<std::sqrt,rowvec>(sum(square(x.each_col() - xnew.col(i)),0)).t();
      }
    }
  }else if(type == "manhattan"){
    for(int i=0;i<nu;++i){
      disa.col(i)=sum(abs(x.each_col() - xnew.col(i)),0).t();
    }
  }
  else stop("Unknown type argument. you have to enter \"euclidean\" or \"manhattan\".");
  return disaa;
}

RcppExport SEXP Rfast_dista(SEXP XnewSEXP,SEXP XSEXP,SEXP sqrSEXP,SEXP typeSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type Xnew(XnewSEXP);
    traits::input_parameter< NumericMatrix >::type X(XSEXP);
    traits::input_parameter< const bool >::type sqr(sqrSEXP);
    traits::input_parameter< const string >::type type(typeSEXP);
    __result = wrap(dista(Xnew,X,sqr,type));
    return __result;

END_RCPP
}

//[[Rcpp::export]]
IntegerMatrix dista_index(NumericMatrix Xnew, NumericMatrix X,const int k,const string type) {
  const int n=X.ncol(),nu=Xnew.ncol(),p=Xnew.nrow();
  mat xnew(Xnew.begin(),p,nu,false),x(X.begin(),p,n,false);
  IntegerMatrix disaa(k,nu);
  imat disa(disaa.begin(),k,nu,false);
  if(type == "euclidean"){
    for(int i=0;i<nu;++i)
      disa.col(i)=get_k_indices(sum(square(x.each_col() - xnew.col(i)),0),k);
  }else if(type == "manhattan"){
    for(int i=0;i<nu;++i)
      disa.col(i)=get_k_indices(sum(abs(x.each_col() - xnew.col(i)),0),k);
  }
  else stop("Unknown type argument. you have to enter \"euclidean\" or \"manhattan\".");
  return disaa;
}

RcppExport SEXP Rfast_dista_index(SEXP XnewSEXP,SEXP XSEXP,SEXP kSEXP,SEXP typeSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type Xnew(XnewSEXP);
    traits::input_parameter< NumericMatrix >::type X(XSEXP);
    traits::input_parameter< const int >::type k(kSEXP);
    traits::input_parameter< const string >::type type(typeSEXP);
    __result = wrap(dista_index(Xnew,X,k,type));
    return __result;

END_RCPP
}

//[[Rcpp::export]]
colvec get_k_values(rowvec x,const int& k){
  sort(x.begin(),x.end());
  return conv_to<colvec>::from(x.subvec(0,k-1));
}

//[[Rcpp::export]]
NumericMatrix dista_values(NumericMatrix Xnew, NumericMatrix X,const int k,const bool sqr,const string type) {
  const int n=X.ncol(),nu=Xnew.ncol();
  mat xnew(Xnew.begin(),Xnew.nrow(),nu,false),x(X.begin(),X.nrow(),n,false);
  NumericMatrix disaa(k,nu);
  mat disa(disaa.begin(),k,nu,false);
  if(type == "euclidean"){
    if(sqr){
      for(int i=0;i<nu;++i)
        disa.col(i)=get_k_values(sum(square(x.each_col() - xnew.col(i)),0),k);
    }else{
      for(int i=0;i<nu;++i){
        disa.col(i)=foreach<std::sqrt,colvec>(get_k_values(sum(square(x.each_col() - xnew.col(i)),0),k));
      }
    }
  }else if(type == "manhattan"){
    for(int i=0;i<nu;++i)
      disa.col(i)=get_k_values(sum(abs(x.each_col() - xnew.col(i)),0),k);
  }
  else stop("Unknown type argument. you have to enter \"euclidean\" or \"manhattan\".");
  return disaa;
}

RcppExport SEXP Rfast_dista_values(SEXP XnewSEXP,SEXP XSEXP,SEXP kSEXP,SEXP sqrSEXP,SEXP typeSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type Xnew(XnewSEXP);
    traits::input_parameter< NumericMatrix >::type X(XSEXP);
    traits::input_parameter< const int >::type k(kSEXP);
    traits::input_parameter< const bool >::type sqr(sqrSEXP);
    traits::input_parameter< const string >::type type(typeSEXP);
    __result = wrap(dista_values(Xnew,X,k,sqr,type));
    return __result;

END_RCPP
}