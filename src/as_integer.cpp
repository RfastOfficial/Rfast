
//Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include "Rfast.h"

using namespace Rcpp;

using std::vector;
using std::string;

//[[Rcpp::export]]
IntegerVector as_integer(SEXP x,const bool sorted,const int init,const bool parallel = false){
  IntegerVector f(LENGTH(x));
  switch(TYPEOF(x)){
    case REALSXP:
      sorted ? as_integer_h_sorted<double>(as< vector<double> >(x),f,init,0,parallel) : as_integer_h<double>(as< vector<double> >(x),f,init,0,parallel);
      break;
    case INTSXP:
      sorted ? as_integer_h_sorted<int>(as< vector<int> >(x),f,init,0,parallel) : as_integer_h<int>(as< vector<int> >(x),f,init,0,parallel);
      break;
    case STRSXP:
      sorted ? as_integer_h_sorted<string>(as< vector<string> >(x),f,init,"",parallel) : as_integer_h<string>(as< vector<string> >(x),f,init,"",parallel);
      break;
    default:
      stop("Wrong type for argument x.\n");
  }
  return f;
}

// // [[Rcpp::export]]
// IntegerVector as_factor(SEXP x){
//   List L;
//   as_integer_h_with_names<double>(as<vector<double>>(x),L,1,0.0);
//   IntegerVector f=L["f"];
//   f.attr("levels") = L["w"];
//   return f;
//   f.attr("class") = "factor";
//   return f;
// }

// // [[Rcpp::export]]
// IntegerVector as_factor2(SEXP x){
//   List L;
//   as_integer_h_with_names<double>(as<vector<double>>(x),L,1,0.0);
//   IntegerVector f=L["f"];
//   f.attr("class") = "factor";
//   f.attr("levels") = as<CharacterVector>(L["w"]);
//   return f;
// }

RcppExport SEXP Rfast_as_integer(SEXP x,SEXP sortedSEXP,SEXP initSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const bool >::type sorted(sortedSEXP);
    traits::input_parameter< const int >::type init(initSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = as_integer(x,sorted,init,parallel);
    return __result;
END_RCPP
}