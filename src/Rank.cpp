
//Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "Rfast.h"
#include "Random.h"

using namespace Rcpp;
using namespace arma;

using std::vector;

static void Rank_random(NumericVector& x,NumericVector& f, const bool parallel = false){
    const unsigned int n=x.size();
    unsigned int i,j=0;
    NumericVector xx=clone(x);
    vector<int> ind(n);

    Random::uniform<Random::integer,true> rng(0,n-1);
		for(i=0;i<n;++i){
			ind[i]=rng();
		}
    int k=0,m,times=0;
    double mn=0.0,v=xx[ind[j]];
    for(i=1;i<n;++i){
        if(v!=xx[ind[i]]){
            times=i-j;
            mn=(j+1+i)*0.5; //mn=mean(seq(j+1,i));
            for(k=j,m=0;m<times;++m)
                f[ind[k++]]=mn;
            j=i;
            v=xx[ind[j]];
        }
    }
    times=i-j;
    mn=(j+i+1)*0.5; //mn=mean(seq(j+1,i));
    for(k=j,m=0;m<times;++m)
        f[ind[k++]]=mn;
}

static void Rank_mean(NumericVector& x,NumericVector& f,const bool descend, const bool parallel = false){
    const int n=x.size();
    int i,j=0;
    NumericVector xx=clone(x);
    vector<int> ind=Order_rank<vector<int>,NumericVector>(xx,descend,false,0,0,parallel);
    int k=0,m,times=0;
    double mn=0.0,v=xx[ind[j]];
    for(i=1;i<n;++i){
        if(v!=xx[ind[i]]){
            times=i-j;
            mn=(j+1+i)*0.5; //mn=mean(seq(j+1,i));
            for(k=j,m=0;m<times;++m)
                f[ind[k++]]=mn;
            j=i;
            v=xx[ind[j]];
        }
    }
    times=i-j;
    mn=(j+i+1)*0.5; //mn=mean(seq(j+1,i));
    for(k=j,m=0;m<times;++m)
        f[ind[k++]]=mn;
}

static void Rank_max(NumericVector& x,NumericVector& f,const bool descend, const bool parallel = false){
    const int n=x.size();
    int i,j=0;
    NumericVector xx=clone(x);
    vector<int> ind=Order_rank<vector<int>,NumericVector>(xx,descend,false,0,0,parallel);
    int k=0,m,times=0;
    double v=xx[ind[0]];
    for(i=1;i<n;++i){
        if(v!=xx[ind[i]]){
            times=i-j;
            for(k=j,m=0;m<times;++m)
                f[ind[k++]]=i;
            j=i;
            v=xx[ind[j]];
        }
    }
    times=i-j;
    for(k=j,m=0;m<times;++m)
        f[ind[k++]]=i;
}

static void Rank_min(NumericVector& x,NumericVector& f,const bool descend, const bool parallel = false){
  const int n=x.size();
  int i,j=0;
  NumericVector xx=clone(x);
  vector<int> ind=Order_rank<vector<int>,NumericVector>(xx,descend,false,0,1,parallel);
  double v=xx[ind[j]];
  f[ind[0]]=1;
  for(i=1;i<n;++i){
    if(v!=xx[ind[i]]){
      j=i;
      v=xx[ind[j]];
    }
    f[ind[i]]=j+1;
  }
}

static void Rank_first(NumericVector& x,NumericVector& f,const bool descend,const bool stable, const bool parallel = false){
  const int n=x.size();
  NumericVector xx=clone(x);
  vector<int> ind=Order_rank<vector<int>,NumericVector>(xx,descend,stable,0,0,parallel);
  for(int i=0;i<n;++i){
    f[ind[i]]=i+1;
  }
}

NumericVector Rank(NumericVector x,string method="average",const bool descend=false,const bool stable=false, const bool parallel = false){
  NumericVector res(x.size());
  if(method == "average"){
    Rank_mean(x,res,descend, parallel);
  }else if(method == "min"){
    Rank_min(x,res,descend, parallel);
  }else if(method == "max"){
    Rank_max(x,res,descend, parallel);
  }else if(method == "first"){
    Rank_first(x,res,descend,stable, parallel);
  }else if(method == "random"){
    Rank_random(x,res, parallel);
  }else{
    stop("Error. Wrong method.");
  }
  return res;
}

RcppExport SEXP Rfast_rank(SEXP xSEXP,SEXP methodSEXP,SEXP descendSEXP,SEXP stableSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< string >::type method(methodSEXP);
    traits::input_parameter< const bool >::type descend(descendSEXP);
    traits::input_parameter< const bool >::type stable(stableSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = Rank(x,method,descend,stable,parallel);
    return __result;
END_RCPP
}
