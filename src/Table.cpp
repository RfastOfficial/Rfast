//Author: Manos Papadakis

//[[Rcpp::plugins(cpp11)]]

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <vector>
#include "Rfast.h"

using namespace Rcpp;

using std::vector;
using std::string;

static void table2_integer(IntegerVector x,IntegerVector y,IntegerMatrix &f){
  const int n=x.size();
  int xx,yy,mn_x,mn_y,mx_x,mx_y;
  min_max<int>(x.begin(),x.end(),mn_x,mx_x);
  min_max<int>(y.begin(),y.end(),mn_y,mx_y);
  f = IntegerMatrix(mx_x-mn_x+1,mx_y-mn_y+1);
  for(int i=0;i<n;++i){
    xx=x[i]-mn_x;
    yy=y[i]-mn_y;
    f(xx,yy)++;
  }
}

static void table2_integer_names(IntegerVector x,IntegerVector y,List &L){
  const int n=x.size();
  int xx,yy,mn_x,mn_y,mx_x,mx_y;
  min_max<int>(x.begin(),x.end(),mn_x,mx_x);
  min_max<int>(y.begin(),y.end(),mn_y,mx_y);
  IntegerMatrix f(mx_x-mn_x+1,mx_y-mn_y+1);
  for(int i=0;i<n;++i){
    xx=x[i]-mn_x;
    yy=y[i]-mn_y;
    f(xx,yy)++;
  }
  L["x"]=seq(mn_x,mx_x);
  L["y"]=seq(mn_y,mx_y);
  L["f"]=f;
}

IntegerMatrix table2_c(SEXP x,SEXP y,const bool rm_zero_col_row=true){
    IntegerMatrix f;
    switch(TYPEOF(x)){
        case REALSXP: 
            table2_like_r<double>(as< vector<double> >(x),as< vector<double> >(y),f,0);
            break;
        case INTSXP: 
            rm_zero_col_row 
            ? table2_like_r<int>(as< vector<int> >(x),as< vector<int> >(y),f,0)
            : table2_integer(IntegerVector(x),IntegerVector(y),f);
            break;
        case STRSXP: 
            table2_like_r<string>(as< vector<string> >(x),as< vector<string> >(y),f,"");
            break;
        default: stop("Wrong type of vector x.");
    }
    return f;
}

List table2_with_names(SEXP x,SEXP y,const bool rm_zero_col_row=true){
    List f;
    switch(TYPEOF(x)){
        case REALSXP: 
            table2_like_r_with_names<double>(as< vector<double> >(x),as< vector<double> >(y),f,0);
            break;
        case INTSXP: 
            rm_zero_col_row 
            ? table2_like_r_with_names<int>(as< vector<int> >(x),as< vector<int> >(y),f,0)
            : table2_integer_names(IntegerVector(x),IntegerVector(y),f);
            break;
        case STRSXP: 
            table2_like_r_with_names<string>(as< vector<string> >(x),as< vector<string> >(y),f,"");
            break;
        default: stop("Wrong type of vector x.");
    }
    return f;
}

vector<int> table_c(SEXP x,const int use_na){
    vector<int> f;
    switch(TYPEOF(x)){
        case REALSXP: 
            f = use_na ? table_use_na<double>(as< vector<double> >(x),use_na) : table_simple<double>(as< vector<double> >(x));
            break;
        case INTSXP: 
            f = use_na ? table_use_na<int>(as< vector<int> >(x),use_na) : table_simple<int>(as< vector<int> >(x));
            break;
        case STRSXP: 
            f = table_simple<string>(as< vector<string> >(x));
            break;
        default: stop("Wrong type of vector x.");
    }
    return f;
}

IntegerVector table_with_names(SEXP x){
  return Rf_isString(x) ? Rcpp::table(as<CharacterVector>(x)) : Rcpp::table(as<NumericVector>(x));
}



RcppExport SEXP Rfast_table_c(SEXP x,SEXP use_naSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type use_na(use_naSEXP);
    __result = table_c(x,use_na);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_table_with_names(SEXP x){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = table_with_names(x);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_table2_c(SEXP x,SEXP y,SEXP rm_zerosSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type rm_zeros(rm_zerosSEXP);
    __result = table2_c(x,y,rm_zeros);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_table2_with_names(SEXP x,SEXP y,SEXP rm_zerosSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type rm_zeros(rm_zerosSEXP);
    __result = table2_with_names(x,y,rm_zeros);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////////////////////////////////


static void table_sign_simple(NumericVector& x,NumericVector& f){
  double v;
  int neg=0,zer=0,pos=0;
  for(NumericVector::iterator xx=x.begin();xx!=x.end();++xx){
    v=*xx;
    if(v>0)
      ++pos;
    else if(v<0)
      ++neg;
    else
      ++zer;
  }
  f[2]=pos;
  f[1]=zer;
  f[0]=neg;
}

static void table_sign_na(NumericVector& x,NumericVector& f){
  double v;
  int neg=0,zer=0,pos=0,nas=0;
  for(NumericVector::iterator xx=x.begin();xx!=x.end();++xx){
    v=*xx;
    if(R_IsNA(v))
      ++nas;
    else if(v>0)
      ++pos;
    else if(v<0)
      ++neg;
    else
      ++zer;
  }
  f[3]=nas;
  f[2]=pos;
  f[1]=zer;
  f[0]=neg;
}

//[[Rcpp::export]]
NumericVector table_sign(NumericVector x,const bool na,const bool names){
  NumericVector f;
  if(na){
    f=NumericVector(4);
    table_sign_na(x,f);
    if(names)
    	f.names()=CharacterVector::create("-1","0","+1","NA");
  }else{
    f=NumericVector(3);
    table_sign_simple(x,f);
    if(names)	
    	f.names()=CharacterVector::create("-1","0","+1");
  }
  return f;
}

RcppExport SEXP Rfast_table_sign(SEXP xSEXP,SEXP naSEXP,SEXP namesSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< const bool >::type na(naSEXP);
    traits::input_parameter< const bool >::type names(namesSEXP);
    __result = table_sign(x,na,names);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////////////