//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include "mn.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

//[[Rcpp::export]]
SEXP cholesky_par(SEXP AA) {
    int i=0,j=0,ni=0,nj=0,n=Rf_ncols(AA);
    SEXP LL=PROTECT(Rf_allocMatrix(REALSXP,n,n));
    double s,*A=REAL(AA),*L=REAL(LL);
    for(j=0;j<n*n;++j)
        L[j]=0;
    for (j = 0; j <n; ++j) {  
        nj=n*j;
        //for (k = 0; k < j; ++k) {s += L[nj + k] * L[nj + k];}
        s = sqrt(A[nj+j] - sum_with<double,square2<double>>(L+nj,L+nj+j));
        L[nj+j] = s;
        s = 1.0/s;
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (i = j+1; i <n; ++i) {
            ni=i*n;
            //for (k = 0; k < j; ++k) {s += L[ni + k] * L[nj + k];}
            L[ni + j] = s * (A[ni + j] - Apply<double,double,mmult<double>,madd<double>>(L+ni,L+ni+j,L+nj));
        }
    }
    UNPROTECT(1);
    return LL;
}

//[[Rcpp::export]]
SEXP cholesky(SEXP AA) {
    int i,j,ni,nj,n=Rf_ncols(AA);
    SEXP LL=PROTECT(Rf_allocMatrix(REALSXP,n,n));
    double s,*A=REAL(AA),*L=REAL(LL);
    for(j=0;j<n*n;++j)
        L[j]=0;
    for (j = 0; j <n; ++j) {  
        nj=n*j;
        //for (k = 0; k < j; ++k) {s += L[nj + k] * L[nj + k];}
        s = sqrt(A[nj+j] - sum_with<double,square2<double>>(L+nj,L+nj+j));
        L[nj+j] = s;
        s = 1.0/s;
        for (i = j+1; i <n; ++i) {
            ni=i*n;
            //for (k = 0; k < j; ++k) {s += L[ni + k] * L[nj + k];}
            L[ni + j] = s * (A[ni + j] - Apply<double,double,mmult<double>,madd<double>>(L+ni,L+ni+j,L+nj));
        }
    }
    UNPROTECT(1);
    return LL;
}

/*
//[[Rcpp::export]]
SEXP cholesky_my(SEXP AA) {
    int i,j,ni,nj,n=Rf_ncols(AA);
    SEXP LL=PROTECT(Rf_allocMatrix(REALSXP,n,n));
    double s,*A=REAL(AA),*L=REAL(LL);
    for(j=0;j<n*n;++j)
        L[j]=0;

    s=sqrt(*A);
    *L=s;
    s=1.0/s;
    for(j=1;j<n;++j){
        nj=n*j;
    	L[nj]=A[nj]*s;
    }

    for(i=1;i<n;++i){
    	ni=n*i;
        //Rcout<<81<<" -> "<<ni<<endl;
    	s=sqrt(A[ni+i]-sum_with<double,square2<double>>(L+ni,L+ni+i));
    	L[ni+i]=s;
    	//Rcout<<84<<" -> "<<ni<<endl;
    	s=1.0/s;
    	for(j=i+1;j<n;++j){
    		nj=n*j;
    	    //Rcout<<88<<" -> "<<nj<<endl;
    		L[nj+i]=s*(A[nj+i] - Apply<double,double,mmult<double>,madd<double>>(L+nj,L+nj+i,L+ni));
    		//Rcout<<81<<endl;
    	}
    }

    UNPROTECT(1);
    return LL;
}


//[[Rcpp::export]]
NumericMatrix cholesky2(NumericMatrix A) {
    int i,j,k,n=A.ncol();
    NumericMatrix L(n,n);
    double s;
    (void)s;
    for (j = 0; j <n; ++j) {  
        s = 0;
        for (k = 0; k < j; ++k) {
          //  s += L(j , k) * L(j , k);
        //    Rcout<<"( "<<j <<" , "<< k<<" ) * ( "<<j <<" , "<< k<<" ) + ";
        }
        //Rcout<<endl;
        //L(j , j) = sqrt(A(j , j) - s);
        //Rcout<<"( "<<j <<" , "<< j<<" )\n";
        for (i = j+1; i <n; ++i) {
            s=0;
            for (k = 0; k < j; ++k) {
              //  s += L(i , k) * L(j , k);                
                Rcout<<"( "<<i <<" , "<< k<<" ) * ( "<<j <<" , "<< k<<" ) + ";
            }
            Rcout<<endl;
            //L(i , j) = (1.0 / L(j , j) * (A(i , j) - s));
            //Rcout<<"( "<<i <<" , "<< j<<" ) ";
        }
        //Rcout<<endl;
    }
    return L;
}

//[[Rcpp::export]]
SEXP init_v(int n){
    SEXP f=Rf_allocMatrix(REALSXP,n,n);
    double *ff=REAL(f);
    for(int i=0;i<n;++i)
        ff[i]=0;
    return f;
}


//[[Rcpp::export]]
SEXP init_v2(int n){
    SEXP f=Rf_allocMatrix(REALSXP,n,n);
    double *ff=REAL(f),*endf=ff+n*n;
    for(;ff!=endf;++ff)
        *ff=0;
    return f;
}*/



RcppExport SEXP Rfast_cholesky(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = wrap(cholesky(x));
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_cholesky_par(SEXP x) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    __result = wrap(cholesky_par(x));
    return __result;
END_RCPP
}
