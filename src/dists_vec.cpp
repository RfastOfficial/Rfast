
//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <string>
#include "mn.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;

using std::string;

static int proper_size(int nrw,int ncl){
  return ncl*(ncl-1)*0.5;
}

//[[Rcpp::export]]
IntegerVector index_dist_vec(const int nrw,const int ncl){
  IntegerVector f(proper_size(nrw,ncl));
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    for(j=i+1;j<ncl;++j,++k){
      f[k]=j+1;
    }
  }
  return f;
}



//[[Rcpp::export]]
List euclidean_dist_vec_ina(NumericMatrix x,const bool sqr){
  const int ncl=x.ncol(),nrw=x.nrow(),n=proper_size(nrw,ncl);
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(n);
  IntegerVector ind_i(n),ind_j(n);
  colvec xv(nrw);
  int i,j,k=0;
  if(sqr){
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j,++k){
        f[k]=sum(square(xx.col(j)-xv));
        ind_i[k]=i+1;
        ind_j[k]=j+1;
      }
    }
  }else{
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j,++k){
        f[k]=std::sqrt(sum(square(xv-xx.col(j))));
        ind_i[k]=i+1;
        ind_j[k]=j+1;
      }
    }
  }
  return List::create(_["i"]=ind_i,_["j"]=ind_j,_["dist"]=f);
}




NumericVector euclidean_dist_vec(NumericMatrix x,const bool sqr){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw);
  int i,j,k=0;
  if(sqr)
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j,++k){
        f[k]=sum(square(xx.col(j)-xv));
      }
    }
  else
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j,++k){
        f[k]=std::sqrt(sum(square(xv-xx.col(j))));
      }
    }
  return f;
}

NumericVector manhattan_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw);
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      f[k]=sum(abs(xv-xx.col(j)));
    }
  }
  return f;
}

NumericVector hellinger_dist_vec(NumericMatrix x,const bool sqr){
  const int ncl=x.ncol(),nrw=x.nrow();
  const double p=1.0/std::sqrt(2.0);
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw);
  int i,j,k=0;
  if(sqr)
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j,++k){
        f[k]=sum(square(xv-xx.col(j)))*0.5;
      }
    }
  else
    for(i=0;i<ncl-1;++i){
      xv=xx.col(i);
      for(j=i+1;j<ncl;++j,++k){
        f[k]=p*std::sqrt(sum(square(xv-xx.col(j))));
      }
    }
  return f;
}

NumericVector max_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw),tmp(nrw);
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      tmp=abs(xv-xx.col(j));
      f[k]=tmp.at(tmp.index_max());
    }
  }
  return f;
}

NumericVector min_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw);
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      xv=abs(xx.col(j)-xv);
      f[k]=xv.at(xv.index_min());
    }
  }
  return f;
}

NumericVector minkowski_dist_vec(NumericMatrix x,const double p){
  const int ncl=x.ncol(),nrw=x.nrow();
  const double p_1=1.0/p;
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw);
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j,++k){      
      f[k]=pow(sum_with<std::pow,colvec>(abs(xv-xx.col(j)),p),p_1);
    }
  }
  return f;
}

NumericVector canberra1_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw),yv(nrw);
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      yv=xx.col(j);
      f[k]=sum(abs((xv-yv)/(xv+yv)));
    }
  }
  return f;
}

NumericVector canberra2_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw),yv(nrw);
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      yv=xx.col(j);
      f[k]=sum(abs(xv-yv)/(abs(xv)-abs(yv)));
    }
  }
  return f;
}



NumericVector total_variation_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw);
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      f[k]=0.5*sum(abs(xv-xx.col(j)));
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericVector kullback_leibler_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  NumericMatrix log_x(nrw,ncl);
  NumericVector f(proper_size(nrw,ncl));
  mat xx(x.begin(),nrw,ncl,false),log_xx(log_x.begin(),nrw,ncl,false);
  colvec xv(nrw),log_xv(nrw);
  int i,j,k=0;
  fill_with<std::log,double*,double*>(x.begin(),x.end(),log_xx.begin());
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    log_xv=log_xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      f[k]=sum((xv-xx.col(j))%(log_xv-log_xx.col(j)));
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericVector bhattacharyya_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  mat xx(x.begin(),nrw,ncl,false);
  NumericVector f(proper_size(nrw,ncl));
  colvec xv(nrw);
  int i,j,k=0;
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      f[k]=sum(sqrt(abs(xv%xx.col(j))));
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericVector itakura_saito_dist_vec(NumericMatrix x){
  const int ncl=x.ncol(),nrw=x.nrow();
  NumericVector f(proper_size(nrw,ncl));
  NumericMatrix log_x(nrw,ncl);
  mat xx(x.begin(),nrw,ncl,false),log_xx(log_x.begin(),nrw,ncl,false);
  colvec xv(nrw),log_xv(nrw);
  int i,j,k=0;
  fill_with<std::log,double*,double*>(x.begin(),x.end(),log_xx.begin());
  for(i=0;i<ncl-1;++i){
    xv=xx.col(i);
    log_xv=log_xx.col(i);
    for(j=i+1;j<ncl;++j,++k){
      f[k]=sum((xv-xx.col(j)-log_xv-log_xx.col(j))-1);
    }
  }
  return f;
}

//[[Rcpp::export]]
NumericMatrix harvesine_dist_vec(NumericMatrix x){
  const int nrw=x.nrow();
  const int nrw_1=nrw-1;
  colvec x0(x.begin(),nrw,false),x1(x.begin()+nrw,nrw,false);
  NumericVector f(proper_size(nrw,nrw));
  mat ff(f.begin(),nrw,nrw,false);
  colvec ind_col(nrw_1);
  colvec a(nrw_1);
  int i;
  
  for(i=0;i<nrw_1;++i){
    span ind(i+1,nrw_1);
    ind_col = x0(ind);
    a = square(sin( 0.5 * (x0[i] -ind_col))) + cos(x0[i]) * (cos(ind_col) % square(sin( 0.5 * (x1[i] - x1(ind)))));
    a = 2 * asin( sqrt(a) );
    ff.insert_rows(a);
  }
  return f;
}


//[[Rcpp::export]]
NumericVector dist_vec(NumericMatrix x,const string method,const bool sqr,const int p){
  if(method == "euclidean" || p==2){
    return euclidean_dist_vec(x,sqr);
  }else if(method == "manhattan" || p==1){
    return manhattan_dist_vec(x);
  }else if(method == "maximum"){
    return max_dist_vec(x);
  }else if(method == "minimum"){
    return min_dist_vec(x);
  }else if(method == "canberra1"){
    return canberra1_dist_vec(x);
  }else if(method == "canberra2"){
    return canberra2_dist_vec(x);
  }else if(method == "minkowski"){
    return minkowski_dist_vec(x,p);
  }else if(method == "bhattacharyya"){
    return bhattacharyya_dist_vec(x);
  }else if(method == "hellinger"){
    return hellinger_dist_vec(x,sqr);
  }else if(method == "total_variation"){
    return total_variation_dist_vec(x);
  }else if(method == "kullback_leibler" || method == "jensen_shannon"){
    return kullback_leibler_dist_vec(x);
  }else if(method == "itakura_saito"){
    return itakura_saito_dist_vec(x);
  }else if(method == "harvesine"){
    return harvesine_dist_vec(x);
  }
  stop("Unsupported Method: %s",method);
}

RcppExport SEXP Rfast_dist_vec(SEXP xSEXP,SEXP methodSEXP,SEXP sqrSEXP,SEXP pSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const string >::type method(methodSEXP);
    traits::input_parameter< const bool >::type sqr(sqrSEXP);
    traits::input_parameter< const int >::type p(pSEXP);
    __result = dist_vec(x,method,sqr,p);
    return __result;
END_RCPP
}