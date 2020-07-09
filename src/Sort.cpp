//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <algorithm>
#include <vector>
#include <string>

using namespace Rcpp;
using std::sort;
using std::vector;
using std::string;
using std::remove_if;

static SEXP Sort_simple(SEXP x,const bool descend){
  SEXP f=PROTECT(Rf_duplicate(x));
  int len=LENGTH(x);
  switch(TYPEOF(x)){
    case INTSXP:{
      int *F=INTEGER(f);
      descend ? sort(F,F+len,std::greater<int>()) : sort(F,F+len);
      break;
    }
    default:{
      double *F=REAL(f);
      descend ? sort(F,F+len,std::greater<double>()) : sort(F,F+len);
      break;
    }
  }
  UNPROTECT(1);
  return f;
}

static NumericVector Sort_na_rm(SEXP x,const bool descend){ // na.rm=NA
  NumericVector f(Rf_duplicate(x));
  const int n=remove_if(f.begin(),f.begin()+f.size(),R_IsNA)-f.begin();
  sort(f.begin(),f.begin()+n);
  return f[Range(0,n-1)];
}

static NumericVector Sort_na_last(SEXP x,const bool descend){ // na.rm=T
  NumericVector f(Rf_duplicate(x));
  const int n=remove_if(f.begin(),f.begin()+f.size(),R_IsNA)-f.begin();
  sort(f.begin(),f.begin()+n);
  for(NumericVector::iterator it=f.begin()+n;it!=f.end();++it)
    *it=NA_REAL;
  return f;
}

SEXP Sort(SEXP x,const bool descend,SEXP na){
  if(Rf_isNull(na))
    return Sort_simple(x,descend);
  else if(R_IsNA(Rf_asReal(na))){
    return Sort_na_rm(x,descend);
  }else{
    return Sort_na_last(x,descend);
  }
  stop("Wrong type of na.last argument.\n");
}

RcppExport SEXP Rfast_Sort(SEXP x,SEXP descendSEXP,SEXP na){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const bool >::type descend(descendSEXP);
    __result = Sort(x,descend,na);
    return __result;
END_RCPP
}


//[[Rcpp::export]]
vector<double> Sort_na_first(vector<double> f,const bool descend){ // na.rm=F
  const int n=remove_if(f.rbegin(),f.rend(),R_IsNA)-f.rbegin();
  descend ? sort(f.end()-n,f.end(),std::greater<double>()) : sort(f.end()-n,f.end());
  for(vector<double>::iterator ff=f.begin();ff!=f.begin()+n;++ff)
    *ff=NA_REAL;
  return f;
}

RcppExport SEXP Rfast_Sort_na_first(SEXP xSEXP,SEXP descendSEXP){
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< vector<double> >::type x(xSEXP);
  traits::input_parameter< const bool >::type descend(descendSEXP);
  __result = Sort_na_first(x,descend);
  return __result;
  END_RCPP
}

//[[Rcpp::export]]
vector<string> sort_string(CharacterVector x,const bool descend){
  vector<string> f(x.begin(),x.end());
  descend ? sort(f.begin(),f.end(),std::greater<string>()) : sort(f.begin(),f.end());
  return f;
}

RcppExport SEXP Rfast_sort_string(SEXP xSEXP,SEXP descendSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< CharacterVector >::type x(xSEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    __result = sort_string(x,descend);
    return __result;
END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////////////////


using std::stable_sort;
using std::greater;

SEXP stable_sort(SEXP x,const bool descend){
  SEXP f=PROTECT(Rf_duplicate(x));
  int len=LENGTH(x);
  switch(TYPEOF(x)){
    case INTSXP:{
      int *F=INTEGER(f);
      descend ? stable_sort(F,F+len,greater<int>()) : stable_sort(F,F+len);
      break;
    }
    default:{
      double *F=REAL(f);
      descend ? stable_sort(F,F+len,greater<double>()) : stable_sort(F,F+len);
      break;
    }
  }
  UNPROTECT(1);
  return f;
}

RcppExport SEXP Rfast_stable_sort(SEXP x,SEXP descendSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const bool >::type descend(descendSEXP);
    __result = stable_sort(x,descend);
    return __result;
END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////////


static void min_max_neg_pos_helper(vector<int> &x,int &mnn,int &mxn,int &mnp,int &mxp,bool &has_pos,bool &has_neg){
    int v;
    for(auto xx=x.begin();xx!=x.end();++xx){
        v=*xx;
        if(v<0){
            has_neg=true;
            if(v<mnn)
                mnn=v;
            else if(v>mxn)
                mxn=v;
        }else{
            has_pos=true;
            if(v>mxp)
                mxp=v;
            else if(v<mnp)
                mnp=v;
        }
    }
}

//[[Rcpp::export]]
vector<int> sort_int(vector<int> x){
    int mnp=INT_MAX,mnn=-1,mxp=0,mxn=INT_MIN;
    bool has_pos=false,has_neg=false;
    min_max_neg_pos_helper(x,mnn,mxn,mnp,mxp,has_pos,has_neg);
    vector<int> pos,f(x.size()),neg;
    vector<int>::iterator a=x.begin(),F=f.begin();
    if(has_pos){
        pos.resize(mxp-mnp+1,0);
    }
    if(has_neg){
        neg.resize(1-mnn+mxn,0);
    }
    if(has_pos && has_neg){
        int aa;
        for(;a!=x.end();++a){
            aa=*a;
            aa<0 ? ++neg[-aa+mxn] : ++pos[aa-mnp];
        }
    }else if(has_pos){
        for(;a!=x.end();++a){
            ++pos[*a-mnp];
        }
    }else{ 
        for(;a!=x.end();++a){
            ++neg[*a-mxn];
        }
    }
    
    
    if(has_neg){
        for(vector<int>::reverse_iterator nr=neg.rbegin();nr!=neg.rend();++nr){
            if(*nr!=0){
                int num=-(neg.rend()-nr-1-mxn),times=*nr;
                for(int i=0;i<times;++i){
                    *F++=num;
                }
            }
        }
    }
    if(has_pos){
        for(a=pos.begin();a!=pos.end();++a){
            if(*a!=0){
                int num=a-pos.begin()+mnp,times=*a;
                for(int i=0;i<times;++i){
                    *F++=num;
                }
            }
        }
    }
    return f;
}

RcppExport SEXP Rfast_sort_int(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< vector<int> >::type x(xSEXP);
    __result = sort_int(x);
    return __result;
END_RCPP
}
