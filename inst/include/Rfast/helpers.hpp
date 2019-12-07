
#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;
using namespace Rcpp;

inline NumericVector col_meds_simple(NumericMatrix& x,const bool parallel=false){
    const int p=x.ncol(),step=x.nrow(),middle=step/2-1;
    NumericVector F(p);
    if(step%2==0){
        if(parallel){
            mat xx(x.begin(),step,p,false);
            colvec ff(F.begin(),p,false),tmpp(step);
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i=0;i<p;++i){
                colvec tmpp=xx.col(i);
                nth_element(tmpp.begin(),tmpp.begin()+middle,tmpp.end());
                ff[i]=(tmpp[middle]+*(min_element(tmpp.begin()+middle+1,tmpp.end())))/2.0;
            }
        }else{
            NumericVector tmp(step);
            for(int i=0;i<p;++i){
                tmp=x.column(i);
                nth_element(tmp.begin(),tmp.begin()+middle,tmp.end());
                F[i]=(tmp[middle]+*(min_element(tmp.begin()+middle+1,tmp.end())))/2.0;
            }
        }
    }else{
        if(parallel){
            mat xx(x.begin(),step,p,false);
            colvec ff(F.begin(),p,false);
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i=0;i<p;++i){
                colvec tmpp=xx.col(i);
                nth_element(tmpp.begin(),tmpp.begin()+middle+1,tmpp.end());
                ff[i]=tmpp[middle+1];
            }
        }else{
            NumericVector tmp(step);
            for(int i=0;i<p;++i){
                tmp=x.column(i);
                nth_element(tmp.begin(),tmp.begin()+middle+1,tmp.end());
                F[i]=tmp[middle+1];
            }
        }
    }
    return F;
}

inline NumericVector col_meds_na_rm(NumericMatrix& x,const bool parallel=false){
    const int p=x.ncol();
    NumericVector F(p);
    if(parallel){
        mat xx(x.begin(),x.nrow(),p,false);
        colvec ff(F.begin(),p,false);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<p;++i){
            colvec tmp=xx.col(i);
            ff[i]=med_helper<colvec>(tmp.begin(),tmp.begin()+(int)(std::remove_if(tmp.begin(),tmp.end(),R_IsNA)-tmp.begin()));
        }
    }else{
        NumericVector tmp(x.nrow());
        for(int i=0;i<p;++i){
            tmp=x.column(i);
            F[i]=med_helper<NumericVector>(tmp.begin(),tmp.begin()+(int)(std::remove_if(tmp.begin(),tmp.end(),R_IsNA)-tmp.begin()));
        } 
    }
    return F;
}


inline rowvec col_meds_simple(mat& x,const bool parallel=false){
    const int p=x.n_cols,step=x.n_rows,middle=step/2-1;
    rowvec F(p);
    if(step%2==0){
        if(parallel){
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i=0;i<p;++i){
                colvec tmp=x.col(i);
                nth_element(tmp.begin(),tmp.begin()+middle,tmp.end());
                F[i]=(tmp[middle]+*(min_element(tmp.begin()+middle+1,tmp.end())))/2.0;
            }
        }else{
            colvec tmp(step);
            for(int i=0;i<p;++i){
                tmp=x.col(i);
                nth_element(tmp.begin(),tmp.begin()+middle,tmp.end());
                F[i]=(tmp[middle]+*(min_element(tmp.begin()+middle+1,tmp.end())))/2.0;
            }
        }
    }else{
        if(parallel){
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i=0;i<p;++i){
                colvec tmp=x.col(i);
                nth_element(tmp.begin(),tmp.begin()+middle+1,tmp.end());
                F[i]=tmp[middle+1];
            }
        }else{
            colvec tmp(step);
            for(int i=0;i<p;++i){
                tmp=x.col(i);
                nth_element(tmp.begin(),tmp.begin()+middle+1,tmp.end());
                F[i]=tmp[middle+1];
            }
        }
    }
    return F;
}

inline rowvec col_meds_na_rm(mat& x,const bool parallel=false){
    const int p=x.n_cols;
    rowvec F(p);
    if(parallel){
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<p;++i){
            colvec tmp=x.col(i);
            F[i]=med_helper<colvec>(tmp.begin(),tmp.begin()+(int)(std::remove_if(tmp.begin(),tmp.end(),R_IsNA)-tmp.begin()));
        }
    }else{
        colvec tmp(x.n_rows);
        for(int i=0;i<p;++i){
            tmp=x.col(i);
            F[i]=med_helper<colvec>(tmp.begin(),tmp.begin()+(int)(std::remove_if(tmp.begin(),tmp.end(),R_IsNA)-tmp.begin()));
        }
    }
    return F;
}

inline NumericVector row_meds_simple(NumericMatrix x,const bool parallel=false){
    const int sz=x.ncol(),p=x.nrow(),middle=sz/2-1;
    NumericVector F(p);
    if(sz%2==0){
        if(parallel){
            mat xx(x.begin(),p,sz,false);
            colvec ff(F.begin(),p,false);
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i=0;i<p;++i){
                rowvec rowii=xx.row(i);
                nth_element(rowii.begin(),rowii.begin()+middle,rowii.end());
                ff[i]=(rowii[middle]+*(min_element(rowii.begin()+middle+1,rowii.end())))/2.0;
            }
        }else{
            NumericVector rowi(sz);
            for(int i=0;i<p;++i){
                rowi=x.row(i);
                nth_element(rowi.begin(),rowi.begin()+middle,rowi.end());
                F[i]=(rowi[middle]+*(min_element(rowi.begin()+middle+1,rowi.end())))/2.0;
            }
        }
    }else{
        if(parallel){
            mat xx(x.begin(),p,sz,false);
            colvec ff(F.begin(),p,false);
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i=0;i<p;++i){
                rowvec rowii=xx.row(i);
                nth_element(rowii.begin(),rowii.begin()+middle,rowii.end());
                ff[i]=rowii[middle+1];
            }
        }else{
            NumericVector rowi(sz);
            for(int i=0;i<p;++i){
                rowi=x.row(i);
                nth_element(rowi.begin(),rowi.begin()+middle,rowi.end());
                F[i]=rowi[middle+1];
            }
        }
    }
    return F;
}

inline colvec row_meds_simple(mat x,const bool parallel=false){
    const int sz=x.n_cols,p=x.n_rows,middle=sz/2-1;
    colvec F(p);
    if(sz%2==0){
        if(parallel){
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i=0;i<p;++i){
                rowvec rowi=x.row(i);
                nth_element(rowi.begin(),rowi.begin()+middle,rowi.end());
                F[i]=(rowi[middle]+*(min_element(rowi.begin()+middle+1,rowi.end())))/2.0;
            }
        }else{
            rowvec rowi(sz);
            for(int i=0;i<p;++i){
                rowi=x.row(i);
                nth_element(rowi.begin(),rowi.begin()+middle,rowi.end());
                F[i]=(rowi[middle]+*(min_element(rowi.begin()+middle+1,rowi.end())))/2.0;
            }
        }
    }else{
        if(parallel){
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i=0;i<p;++i){
                rowvec rowi=x.row(i);
                nth_element(rowi.begin(),rowi.begin()+middle,rowi.end());
                F[i]=rowi[middle+1];
            }
        }else{
            rowvec rowi(sz);
            for(int i=0;i<p;++i){
                rowi=x.row(i);
                nth_element(rowi.begin(),rowi.begin()+middle,rowi.end());
                F[i]=rowi[middle+1];
            }
        }
    }
    return F;
}

inline colvec row_meds_na_rm(mat& x,const bool parallel=false){
    const int p=x.n_rows;
    colvec F(p);
    if(parallel){
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<p;++i){
            rowvec tmp=x.row(i);
            F[i]=med_helper<rowvec>(tmp.begin(),tmp.begin()+(int)(std::remove_if(tmp.begin(),tmp.end(),R_IsNA)-tmp.begin()));
        }
    }else{
        rowvec tmp(x.n_cols);
        for(int i=0;i<p;++i){
            tmp=x.row(i);
            F[i]=med_helper<rowvec>(tmp.begin(),tmp.begin()+(int)(std::remove_if(tmp.begin(),tmp.end(),R_IsNA)-tmp.begin()));
        }
    }
    return F;
}

inline NumericVector row_meds_na_rm(NumericMatrix& x,const bool parallel=false){
    const int p=x.nrow();
    NumericVector F(p);
    if(parallel){
        mat xx(x.begin(),x.nrow(),p,false);
        colvec ff(F.begin(),p,false);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<p;++i){
            rowvec tmp=xx.row(i);
            ff[i]=med_helper<rowvec>(tmp.begin(),tmp.begin()+(int)(std::remove_if(tmp.begin(),tmp.end(),R_IsNA)-tmp.begin()));
        }
    }else{
        NumericVector tmp(x.ncol());
        for(int i=0;i<p;++i){
            tmp=x.row(i);
            F[i]=med_helper<rowvec>(tmp.begin(),tmp.begin()+(int)(std::remove_if(tmp.begin(),tmp.end(),R_IsNA)-tmp.begin()));
        }
    }
    return F;
}

inline long long int get_current_nanoseconds(){
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

#endif
