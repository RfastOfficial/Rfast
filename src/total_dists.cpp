
//Author: Manos Papadakis

//[[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;


double total_euclidean_dist(NumericMatrix x,const bool sqr){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw);
  double a=0;
  int i,j;
  if(sqr){
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j){
        a+=sum(square(xx.col(j)-xv));
      }
    }
  }else{
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j){
        a+=std::sqrt(sum(square(xv-xx.col(j))));
      }
    }
  }
  return a;
}

double total_manhattan_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw);
  double a=0;
  int i,j;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j){
      a+=sum(abs(xv-xx.col(j)));
    }
  }
  return a;
}

double total_hellinger_dist(NumericMatrix x,const bool sqr){
  const int ncl=x.ncol(),nrw=x.nrow();
  const double p=1.0/std::sqrt(2.0);
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw);
  double a=0;
  int i,j;
  if(sqr){
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j){
        a+=sum(square(xv-xx.col(j)))*0.5;
      }
    }
  }else{
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j){
        a+=p*std::sqrt(sum(square(xv-xx.col(j))));
      }
    }
  }
  return a;
}

double total_max_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw),tmp(nrw);
  double a=0;
  int i,j;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j){
      tmp=abs(xv-xx.col(j));
      a+=tmp.at(tmp.index_max());
      
    }
  }
  return a;
}

double total_min_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw);
  double a=0;
  int i,j;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j){
      xv=abs(xx.col(j)-xv);
      a+=xv.at(xv.index_min());
      
    }
  }
  return a;
}

double total_minkowski_dist(NumericMatrix x,const double p){
  const int ncl=x.ncol(),nrw=x.nrow();
  const double p_1=1.0/p;
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw);
  double a=0;
  int i,j;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j){      
      a+=pow(sum_with<std::pow,colvec>(abs(xv-xx.col(j)),p),p_1);
    }
  }
  return a;
}

double total_canberra1_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw),yv(nrw);
  double a=0;
  int i,j;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j){
      yv=xx.col(j);
      a+=sum(abs((xv-yv)/(xv+yv)));
    }
  }
  return a;
}

double total_canberra2_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw),yv(nrw);
  double a=0;
  int i,j;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j){
      yv=xx.col(j);
      a+=sum(abs(xv-yv)/(abs(xv)-abs(yv)));
      
    }
  }
  return a;
}



double total_total_variation_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw);
  double a=0;
  int i,j;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j){
      a+=0.5*sum(abs(xv-xx.col(j)));
    }
  }
  return a;
}


double total_kullback_leibler_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  NumericMatrix log_x(nrw,ncl);
  mat xx(x.begin(),nrw,ncl,false),log_xx(log_x.begin(),nrw,ncl,false);
  colvec xv(nrw),log_xv(nrw);
  double a=0;
  int i,j;
  fill_with<std::log,double*,double*>(x.begin(),x.end(),log_xx.begin());
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    log_xv=log_xx.col(i);
    for(j=i+1;j<ncl;++j){
      a+=sum((xv-xx.col(j))%(log_xv-log_xx.col(j)));
    }
  }
  return a;
}


double total_bhattacharyya_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  colvec xv(nrw);
  double a=0;
  int i,j;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j){
      a+=sum(sqrt(abs(xv%xx.col(j))));
    }
  }
  return a;
}


double total_itakura_saito_dist(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  NumericMatrix log_x(nrw,ncl);
  mat xx(x.begin(),nrw,ncl,false),log_xx(log_x.begin(),nrw,ncl,false);
  colvec xv(nrw),log_xv(nrw);
  double a=0;
  int i,j;
  fill_with<std::log,double*,double*>(x.begin(),x.end(),log_xx.begin());
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    log_xv=log_xx.col(i);
    for(j=i+1;j<ncl;++j){
      a+=sum((xv-xx.col(j)-log_xv-log_xx.col(j))-1);
    }
  }
  return a;
}


double total_dists(NumericMatrix x,const string method,const bool sqr,const int p){
  if(method == "euclidean" || p==2){
    return total_euclidean_dist(x,sqr);
  }else if(method == "manhattan" || p==1){
    return total_manhattan_dist(x);
  }else if(method == "maximum"){
    return total_max_dist(x);
  }else if(method == "minimum"){
    return total_min_dist(x);
  }else if(method == "canberra1"){
    return total_canberra1_dist(x);
  }else if(method == "canberra2"){
    return total_canberra2_dist(x);
  }else if(method == "minkowski"){
    return total_minkowski_dist(x,p);
  }else if(method == "bhattacharyya"){
    return total_bhattacharyya_dist(x);
  }else if(method == "hellinger"){
    return total_hellinger_dist(x,sqr);
  }else if(method == "total_variation"){
    return total_total_variation_dist(x);
  }else if(method == "kullback_leibler" || method == "jensen_shannon"){
    return total_kullback_leibler_dist(x);
  }else if(method == "itakura_saito"){
    return total_itakura_saito_dist(x);
  }
  stop("Unsupported Method: %s",method);
}

RcppExport SEXP Rfast_total_dists(SEXP xSEXP,SEXP methodSEXP,SEXP sqrSEXP,SEXP pSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const string >::type method(methodSEXP);
    traits::input_parameter< const bool >::type sqr(sqrSEXP);
    traits::input_parameter< const int >::type p(pSEXP);
    __result = total_dists(x,method,sqr,p);
    return __result;
END_RCPP
}