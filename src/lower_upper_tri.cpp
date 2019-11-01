//Author: Manos Papadakis


#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector upper_tri(NumericMatrix x,const bool dg){
  const int ncl=x.ncol(),nrw=x.nrow();
  int i,j;
  NumericVector f = (ncl<nrw) ? NumericVector(ncl*(nrw-1)*0.5 + (dg ? nrw : 0)) : 
                                NumericVector(nrw*(nrw-1)*0.5 + (dg ? nrw : 0));
  NumericVector::iterator ff=f.begin();
  
  if(dg){
    for(i=0;i<ncl;++i){
      for(j=0;j<=i;++j,++ff){
          *ff=x(j,i);
      }
    }
  }else{
    for(i=1;i<ncl;++i){
      for(j=0;j<i;++j,++ff){       
          *ff=x(j,i);
      }           
    }
  }
  return f;
}


double sum_upper_tri(NumericMatrix x,const bool dg){
  const int ncl=x.ncol();
  int i,j;
  double s=0.0;
  if(dg){
    for(i=0;i<ncl;++i){
      for(j=0;j<=i;++j){
          s+=x(j,i);
      }
    }
  }else{
    for(i=1;i<ncl;++i){
      for(j=0;j<i;++j){
          s+=x(j,i);
      }
    }
  }
  return s;
}


//[[Rcpp::export]]
LogicalMatrix upper_tri_b(int nrw, int ncl,const bool dg){
  int i,j;
  LogicalMatrix f(nrw,ncl);
  if(dg){
    for(i=0;i<ncl;++i){
        for(j=0;j<=i;++j){
          f(j,i)=true;
        }
    }
  }else{
    for(i=1;i<ncl;++i){
        for(j=0;j<i;++j){
          f(j,i)=true;
        }
    }
  }
  return f;
}

RcppExport SEXP Rfast_upper_tri(SEXP xSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type dg(dgSEXP);
    __result = wrap(upper_tri(x,dg));
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_upper_tri_b(SEXP nclSEXP, SEXP nrwSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type ncl(nclSEXP);
    traits::input_parameter< const int >::type nrw(nrwSEXP);
    traits::input_parameter< const bool >::type dg(dgSEXP);
    __result = wrap(upper_tri_b(nrw,ncl,dg));
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_sum_upper_tri(SEXP xSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type dg(dgSEXP);
    __result = wrap(sum_upper_tri(x,dg));
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////////////////////


NumericVector lower_tri(NumericMatrix x,const bool dg){
  const int ncl=x.ncol(),nrw=x.nrow();
  int i,j;
  NumericVector f = (ncl<nrw) ? NumericVector(ncl*(nrw-1)*0.5 + (dg ? nrw : 0)) : 
                                NumericVector(nrw*(nrw-1)*0.5 + (dg ? nrw : 0));
  NumericVector::iterator ff=f.begin();
  if(dg){
    for(i=0;i<ncl;++i){
        for(j=i;j<nrw;++j,++ff){
          *ff=x(j,i);
        }
      }
  }else{
    for(i=0;i<ncl;++i){
        for(j=i+1;j<nrw;++j,++ff){
          *ff=x(j,i);
        }
    }
  }
  return f;
}


double sum_lower_tri(NumericMatrix x,const bool dg){
  const int ncl=x.ncol(),nrw=x.nrow();
  int i,j;
  double s=0.0;
  if(dg){
    for(i=0;i<ncl;++i){
      for(j=i;j<nrw;++j){
        s+=x(j,i);
      }
    }
  }else{
    for(i=0;i<ncl;++i){
      for(j=i+1;j<nrw;++j){
        s+=x(j,i);
      }
    }
  }
  return s;
}

LogicalMatrix lower_tri_b(const int nrw,const int ncl,const bool dg){
  int i,j;
  LogicalMatrix f(nrw,ncl);
  if(dg){
    for(i=0;i<ncl;++i){
      for(j=i;j<nrw;++j){
        f(j,i)=true;
      }
    }
  }else{
    for(i=0;i<ncl;++i){
      for(j=i+1;j<nrw;++j){
        f(j,i)=true;
      }
    }
  }
  return f;
}

RcppExport SEXP Rfast_lower_tri(SEXP xSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type dg(dgSEXP);
    __result = wrap(lower_tri(x,dg));
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_lower_tri_b(SEXP nclSEXP, SEXP nrwSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type ncl(nclSEXP);
    traits::input_parameter< const int >::type nrw(nrwSEXP);
    traits::input_parameter< const bool >::type dg(dgSEXP);
    __result = wrap(lower_tri_b(nrw,ncl,dg));
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_sum_lower_tri(SEXP xSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type dg(dgSEXP);
    __result = wrap(sum_lower_tri(x,dg));
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////////////////////


NumericVector lower_tri_assign(NumericMatrix x,NumericVector v,const bool dg=false){
    const int ncl=x.ncol(),nrw=x.nrow();
    int i,j;
    NumericVector f = clone(x);
    NumericVector::iterator vv=v.begin();
    if(dg){
        for(i=0;i<ncl;++i){
            for(j=i;j<nrw;++j){
                f(j,i)=*vv++;
            }
        }
    }else{
        for(i=0;i<ncl;++i){
            for(j=i+1;j<nrw;++j){
                f(j,i)=*vv++;
            }
        }
    }
    return f;
}

RcppExport SEXP Rfast_lower_tri_assign(SEXP xSEXP,SEXP vSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type v(vSEXP);
    traits::input_parameter< const double >::type dg(dgSEXP);
    __result = wrap(lower_tri_assign(x,v,dg));
    return __result;
END_RCPP
}

NumericVector upper_tri_assign(NumericMatrix x,NumericVector v,const bool dg){
  const int ncl=x.ncol();
  int i,j;
  NumericVector f = clone(x);
  NumericVector::iterator vv=v.begin();  
  if(dg){
    for(i=0;i<ncl;++i){
      for(j=0;j<=i;++j){
          f(j,i)=*vv++;
      }
    }
  }else{
    for(i=1;i<ncl;++i){
      for(j=0;j<i;++j){       
          f(j,i)=*vv++;
      }           
    }
  }
  return f;
}

RcppExport SEXP Rfast_upper_tri_assign(SEXP xSEXP,SEXP vSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type v(vSEXP);
    traits::input_parameter< const double >::type dg(dgSEXP);
    __result = wrap(upper_tri_assign(x,v,dg));
    return __result;
END_RCPP
}
