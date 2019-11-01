//Author: Manos Papadakis


#include <RcppArmadillo.h>
#include "Rfast.h"

using namespace Rcpp;

using std::string;

//[[Rcpp::export]]
SEXP apply_condition(SEXP x,string method,string oper,int val){
    const int p=Rf_ncols(x),n=Rf_nrows(x);
    SEXP f=Rf_allocVector(INTSXP,p);
    int *ff=INTEGER(f),*xx=INTEGER(x),*endx=xx+LENGTH(x);
    if(method == "+"){
        if(oper == ">"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater<int>,madd<int>>(xx,xx+n,val);
            }
        }else if(oper == ">="){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater_eq<int>,madd<int>>(xx,xx+n,val);
            }
        }else if(oper == "<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless<int>,madd<int>>(xx,xx+n,val);
            }
        }else if(oper == "=<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless_eq<int>,madd<int>>(xx,xx+n,val);
            }
        }else{
        stop("Error: Unsupported operation.");
    }
    }else if(method == "-"){
        if(oper == ">"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater<int>,mdiff<int>>(xx,xx+n,val);
            }
        }else if(oper == ">="){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater_eq<int>,mdiff<int>>(xx,xx+n,val);
            }
        }else if(oper == "<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless<int>,mdiff<int>>(xx,xx+n,val);
            }
        }else if(oper == "=<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless_eq<int>,mdiff<int>>(xx,xx+n,val);
            }
        }else{
        stop("Error: Unsupported operation.");
    }
    }else if(method == "min"){
        if(oper == ">"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater<int>,mmin<int>>(xx,xx+n,val);
            }
        }else if(oper == ">="){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater_eq<int>,mmin<int>>(xx,xx+n,val);
            }
        }else if(oper == "<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless<int>,mmin<int>>(xx,xx+n,val);
            }
        }else if(oper == "=<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless_eq<int>,mmin<int>>(xx,xx+n,val);
            }
        }else{
        stop("Error: Unsupported operation.");
    }
    }else if(method == "max"){
        if(oper == ">"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater<int>,mmax<int>>(xx,xx+n,val);
            }
        }else if(oper == ">="){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater_eq<int>,mmax<int>>(xx,xx+n,val);
            }
        }else if(oper == "<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless<int>,mmax<int>>(xx,xx+n,val);
            }
        }else if(oper == "=<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless_eq<int>,mmax<int>>(xx,xx+n,val);
            }
        }else{
        stop("Error: Unsupported operation.");
    }
    }else if(method == "*"){
        if(oper == ">"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater<int>,mmult<int>>(xx,xx+n,val);
            }
        }else if(oper == ">="){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater_eq<int>,mmult<int>>(xx,xx+n,val);
            }
        }else if(oper == "<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless<int>,mmult<int>>(xx,xx+n,val);
            }
        }else if(oper == "=<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless_eq<int>,mmult<int>>(xx,xx+n,val);
            }
        }else{
        stop("Error: Unsupported operation.");
    }
    }else if(method == "/"){
        if(oper == ">"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater<int>,mdiv<int>>(xx,xx+n,val);
            }
        }else if(oper == ">="){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mgreater_eq<int>,mdiv<int>>(xx,xx+n,val);
            }
        }else if(oper == "<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless<int>,mdiv<int>>(xx,xx+n,val);
            }
        }else if(oper == "=<"){
            for(;xx!=endx;xx+=n,++ff){
                *ff=Apply_helper<mless_eq<int>,mdiv<int>>(xx,xx+n,val);
            }
        }else{
        stop("Error: Unsupported operation.");
    }
    }else{
        stop("Error: Unsupported method.");
    }
    return f;
}


RcppExport SEXP Rfast_apply_condition(SEXP x,SEXP methodSEXP,SEXP operSEXP,SEXP valSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< string >::type method(methodSEXP);
    traits::input_parameter< string >::type oper(operSEXP);
    traits::input_parameter< int >::type val(valSEXP);
    __result = wrap(apply_condition(x,method,oper,val));
    return __result;
END_RCPP
}
