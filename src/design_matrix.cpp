//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"

using namespace Rcpp;

IntegerMatrix design_matrix(CharacterVector x,bool ones_c) {
  int i=0;
  const int n=x.size();
  CharacterVector tmp=sort_unique(x);
  CharacterVector::iterator xx=x.begin(),leksi_bg,leksi_en;
  IntegerMatrix Final(n,tmp.size());
  for(leksi_bg=tmp.begin(),leksi_en=tmp.end(),i=0;xx!=x.end();++xx,++i)
    Final(i,lower_bound(leksi_bg,leksi_en,*xx)-leksi_bg)=1;
  if(ones_c){
    IntegerVector ones(n,true);
    Final.column(0)=ones;
  }
  return Final;
}

umat design_matrix_big(DataFrame x) {
  unsigned int i,n=x.length(),last=1,j,nrw;
  umat dm,F;
  dm=design_matrix_helper_big(x(0));
  nrw=dm.n_rows;
  F.reshape(nrw,dm.n_cols);
  for(j=1;j<dm.n_cols;++j)
    F.col(last++)=dm.col(j);
  for(i=1;i<n;++i){
    dm=design_matrix_helper_big(x(i));
    last=F.n_cols;
    F.reshape(nrw,last+dm.n_cols-1); //-1 giati afairo mia stili
    for(j=1;j<dm.n_cols;++j)
      F.col(last++)=dm.col(j);
  }
  F.col(0).fill(1);
  return F;
}

//the model.matrix form R for many columns
RcppExport SEXP Rfast_design_matrix_big(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< DataFrame >::type x(xSEXP);
    __result = design_matrix_big(x);
    return __result;
END_RCPP
}

//the model.matrix form R but by collumn
RcppExport SEXP Rfast_design_matrix(SEXP xSEXP,SEXP onesSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< CharacterVector >::type x(xSEXP);
    traits::input_parameter< bool >::type ones(onesSEXP);
    __result = design_matrix(x,ones);
    return __result;
END_RCPP
}