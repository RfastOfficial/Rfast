
//Author: Manos Papadakis

#include <RcppArmadillo.h>
using namespace Rcpp;

static bool check_if_d(double& x){
  x=std::abs(x);
  return (x-int(x))!=0;
}

bool check_all_ints(NumericVector x){
  bool s=true;
  for(double& v : x)
    if(check_if_d(v)){
      s=false;
      break;
    }
  return s;
}

RcppExport SEXP Rfast_is_integer(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = check_all_ints(x);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////


using std::basic_string;

bool is_element(NumericVector x,double el){
  NumericVector::iterator a=x.begin();
  for(;a!=x.end() && *a!=el;++a);
  return *a==el;
}

RcppExport SEXP Rfast_is_element(SEXP xSEXP,SEXP elSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< double >::type el(elSEXP);
    __result = is_element(x,el);
    return __result;
END_RCPP
}

//[[Rcpp::export]]
bool is_element_string(CharacterVector x,basic_string<char> el){
  CharacterVector::iterator a=x.begin();
  for(;a!=x.end() && *a!=el;++a);
  return *a==el;
}

RcppExport SEXP Rfast_is_element_string(SEXP xSEXP,SEXP elSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< CharacterVector >::type x(xSEXP);
    traits::input_parameter< basic_string<char> >::type el(elSEXP);
    __result = is_element_string(x,el);
    return __result;
END_RCPP
}


/////////////////////////////////////////////