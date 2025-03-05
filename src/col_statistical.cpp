//Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>

using namespace Rcpp;
using namespace arma;

SEXP col_diffs(SEXP x){
  const int n=Rf_nrows(x),p=Rf_ncols(x);
  SEXP f=Rf_allocMatrix(REALSXP,n,p-1);
  double *ff=REAL(f),*xx=REAL(x),*l=xx+n,*end=ff+Rf_length(f);
  for(;ff!=end;++ff,++l,++xx)
    *ff=*l-*xx;
  return f;
}

RcppExport SEXP Rfast_col_diffs(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = col_diffs(x);
    return __result;
END_RCPP
}


///////////////////////////////////////////////////////////////////////////


//[[Rcpp::plugins(cpp11)]]
template<class T>
static double group_sum_tabulate_div(colvec& x,T *kk,T mn,T mx){
  colvec val_m(mx-mn+1,fill::zeros);
  colvec val_ni(mx-mn+1,fill::zeros);
  colvec::iterator xx=x.begin();
  int index;
  double s=0;
  for(;xx!=x.end();++xx,++kk){
  	index=*kk-mn;
    val_m[index]+=*xx;
    val_ni[index]++;
  }
  for(auto mm=val_m.begin(),nni=val_ni.begin();mm!=val_m.end();++mm,++nni){
    if(*mm!=0){
      s+=(*mm * *mm)/ *nni;
    }
  }
  return s;
}



//[[Rcpp::export]]
NumericVector col_anovas(NumericVector Y,IntegerMatrix X) {
  const int nrw=X.nrow(),ncl=X.ncol();
  NumericVector a(ncl),k(ncl);
  colvec kk(k.begin(),ncl,false),y(Y.begin(),nrw,false);
  imat x(reinterpret_cast<imat::elem_type*>(X.begin()),nrw,ncl,false);
  irowvec mx=max(x,0),mn=min(x,0);
  imat::elem_type *startx=x.begin(),*endx=startx+x.n_elem,*tend=startx+nrw,*mxx=mx.begin(),*mnn=mn.begin();
  double *aa=a.begin();
  for (;startx!=endx;++aa,++mxx,++mnn) {
    *aa=group_sum_tabulate_div<imat::elem_type>(y,startx,*mnn,*mxx);
    startx=tend;
    tend+=nrw;
  }
  return a;
}


RcppExport SEXP Rfast_col_anovas(SEXP ySEXP,SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    __result = col_anovas(y,x);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////////////////