// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <algorithm>
#include "mn.h"


using namespace Rcpp;

using std::vector;
using std::string;
using std::sort;

List hash2list(List x,bool sorting){
	CharacterVector nam=x.names();
	vector<string> lnames(nam.begin(),nam.end());
	vector<string>::iterator names=lnames.begin();
	NumericVector val,rr;
	int sz=lnames.size();
	List f(sz);
	List::iterator ff=f.begin();
	if(sorting){
		for(;ff!=f.end();++ff,++names){
			rr=toNumbers(*names," ");
			sort(rr.begin(),rr.end());
			val=x[*names];
			rr.push_back(val[0]);
			*ff=rr;
		}
	}else{
	  	for(;ff!=f.end();++ff,++names){
			rr=toNumbers(*names," ");
			val=x[*names];
			rr.push_back(val[0]);
			*ff=rr;
	  	}
	}
  return f;
}

RcppExport SEXP Rfast_hash2list(SEXP xSEXP,SEXP sortingSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< List >::type x(xSEXP);
    traits::input_parameter< bool >::type sorting(sortingSEXP);
    __result = hash2list(x,sorting);
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////////////////////


using std::vector;
using std::string;
using std::find;

double hash_find(List x,string value){
  vector<string> nam=x.names();
  if(find(nam.begin(),nam.end(),value)!=nam.end())
    return as<double>(x[value]);
  return 0.0;
}

RcppExport SEXP Rfast_hash_find(SEXP xSEXP,SEXP valueSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< List >::type x(xSEXP);
    traits::input_parameter< string >::type value(valueSEXP);
    __result = hash_find(x,value);
    return __result;
END_RCPP
}

///////////////////////////////////////////////////////////////////////////


using std::string;

List Hash_list(CharacterVector key,NumericVector x){
  CharacterVector::iterator k=key.begin();
  NumericVector::iterator xx=x.begin();
  List H_M;
  string m;
  for(;k!=key.end();++k,++xx){
    m=as<string>(*k);
    H_M[m]=*xx;
  }
  return H_M;
}

RcppExport SEXP Rfast_Hash_list(SEXP keySEXP,SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< CharacterVector >::type key(keySEXP);
    traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = Hash_list(key,x);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////////////////////



static bool if_exists_inside(string x,string& value,const string& spliter){
    x+=spliter;
    const char *split=spliter.c_str();
    char *token = std::strtok(&x[0], split);
    while (token != nullptr) {
        if(string(token) == value){
            return true;
        }
        token = std::strtok(nullptr, split);
    }
    return false;
}


string Hash_key_multi(Environment x,string value,const string sep = " "){
    CharacterVector keys=x.ls(true);
    string res;
    for(auto& key : keys){
        res = key;
        if(if_exists_inside(res,value,sep)){
            return res;
        }
    }
    return "";
}

RcppExport SEXP Rfast_Hash_key_multi(SEXP xSEXP,SEXP valueSEXP,SEXP sepSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< Environment >::type x(xSEXP);
    traits::input_parameter< string >::type value(valueSEXP);
    traits::input_parameter< const string >::type sep(sepSEXP);
    __result = Hash_key_multi(x,value,sep);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////////////////////
