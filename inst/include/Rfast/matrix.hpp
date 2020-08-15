
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <thread>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "templates.h"
#include "helpers.hpp"
#include "vector.hpp"

namespace Rfast {
	namespace matrix {

		using namespace Rcpp;
		using namespace arma;
		using namespace std;
		using namespace chrono;

		inline NumericMatrix transpose(NumericMatrix x){
			const int p=x.ncol(),n=x.nrow();
			NumericMatrix f = p==n ? clone(x) : NumericMatrix(p,n);
			if(p==n){
				for(int i=1;i<p;++i){
					for(int u=0;u<i;++u){
						swap(f(u,i),f(i,u));
					}
				}
			}else{
				mat ff(f.begin(),p,n,false),xx(x.begin(),n,p,false);
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for(int i=0;i<p;++i){
					ff.row(i)=xx.col(i).t();
				}
			}
			return f;
		}

		inline mat transpose(mat x){
			const int p=x.n_cols,n=x.n_rows;
			mat f;
			if(p==n){
				f=x;
				for(int i=1;i<p;++i){
					for(int u=0;u<i;++u){
						swap(f(u,i),f(i,u));
					}
				}
			}else{
				f=mat(p,n);
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for(int i=0;i<p;++i){
					f.row(i)=x.col(i).t();
				}
			}
			return f;
		}

		inline NumericMatrix matrix_multiplication(NumericMatrix X,NumericMatrix Y,const bool tx=false,const bool ty=false){
			int p,n;

			if(!tx){
				n=X.nrow();
				p= ty ? Y.nrow() : Y.ncol();
			}else if(tx){
				n=X.ncol();
				p=Y.ncol();
			}
			NumericMatrix C(n,p);
			mat CC(C.begin(),n,p,false),x(X.begin(),X.nrow(),X.ncol(),false),y(Y.begin(),Y.nrow(),Y.ncol(),false);
			colvec yi;

			if(!tx and !ty){ // matrix multiplication
				mat xx=transpose(x);
				for(int i=0;i<p;++i){
				  yi=y.col(i);
				  #ifdef _OPENMP
				  #pragma omp parallel for
				  #endif
				  for(int j=0;j<n;++j){
				    CC(j,i)=dot(xx.col(j),yi);
				  }
				}
			}else if(tx){ // crossprod
				for(int i=0;i<p;++i){
				  yi=y.col(i);
				  #ifdef _OPENMP
				  #pragma omp parallel for
				  #endif
				  for(int j=0;j<n;++j){
				    CC(j,i)=dot(x.col(j),yi);
				  }
				}
			}else{ // tcrossprod
				mat yy=transpose(y);
				mat xx=transpose(x);
				for(int i=0;i<p;++i){
				  yi=yy.col(i);
				  #ifdef _OPENMP
				  #pragma omp parallel for
				  #endif
				  for(int j=0;j<n;++j){
				    CC(j,i)=dot(xx.col(j),yi);
				  }
				}
			}
			return C;
		}

		inline mat matrix_multiplication(mat x,mat y,const bool tx=false,const bool ty=false){
			int p,n;

			if(!tx){
				n=x.n_rows;
				p= ty ? x.n_rows : x.n_cols;
			}else if(tx){
				n=x.n_cols;
				p=x.n_cols;
			}
			mat C(n,p);
			colvec yi;

			if(!tx and !ty){ // matrix multiplication
				mat xx=Rfast::matrix::transpose(x);
				for(int i=0;i<p;++i){
					yi=y.col(i);
					#ifdef _OPENMP
					#pragma omp parallel for
					#endif
					for(int j=0;j<n;++j){
					C(j,i)=dot(xx.col(j),yi);
					}
				}
			}else if(tx){ // crossprod
				for(int i=0;i<p;++i){
				  yi=y.col(i);
				  #ifdef _OPENMP
				  #pragma omp parallel for
				  #endif
				  for(int j=0;j<n;++j){
				    C(j,i)=dot(x.col(j),yi);
				  }
				}
			}else{ // tcrossprod
				mat yy=Rfast::matrix::transpose(y);
				mat xx=Rfast::matrix::transpose(x);
				for(int i=0;i<p;++i){
				  yi=yy.col(i);
				  #ifdef _OPENMP
				  #pragma omp parallel for
				  #endif
				  for(int j=0;j<n;++j){
				    C(j,i)=dot(xx.col(j),yi);
				  }
				}
			}
			return C;
		}

		inline NumericMatrix colSort(NumericMatrix x,const bool descend=false,const bool stable=false,const bool parallel=false){
			const int n=x.nrow(),p=x.ncol();
			NumericMatrix f(n,p);
			mat xx(x.begin(),n,p,false),ff(f.begin(),n,p,false);
			if(descend){
				if(stable){
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<p;++i){
							colvec coli=xx.col(i);
							stable_sort(coli.begin(),coli.end(),greater<double>());
							ff.col(i)=coli;
					    }
					}else{
						colvec coli(n);
						for(int i=0;i<p;++i){
							coli=xx.col(i);
							stable_sort(coli.begin(),coli.end(),greater<double>());
							ff.col(i)=coli;
					    }
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<p;++i){
							colvec coli=xx.col(i);
							sort(coli.begin(),coli.end(),greater<double>());
							ff.col(i)=coli;
					    }
					}else{
						colvec coli(n);
						for(int i=0;i<p;++i){
							coli=xx.col(i);
							sort(coli.begin(),coli.end(),greater<double>());
							ff.col(i)=coli;
					    }
					}
				}
			}
			else{
				if(stable){
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<p;++i){
							colvec coli=xx.col(i);
							stable_sort(coli.begin(),coli.end());
							ff.col(i)=coli;
					    }
					}else{
						colvec coli(n);
						for(int i=0;i<p;++i){
							coli=xx.col(i);
							stable_sort(coli.begin(),coli.end());
							ff.col(i)=coli;
					    }
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<p;++i){
							colvec coli=xx.col(i);
							sort(coli.begin(),coli.end());
							ff.col(i)=coli;
					    }
					}else{
						colvec coli(n);
						for(int i=0;i<p;++i){
							coli=xx.col(i);
							sort(coli.begin(),coli.end());
							ff.col(i)=coli;
					    }
					}
				}
			}  
			return f;
		}

		inline mat colSort(mat x,const bool descend=false,const bool stable=false,const bool parallel=false){
			const int n=x.n_rows,p=x.n_cols;
			mat f(n,p);
			if(descend){
				if(stable){
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<p;++i){
							colvec coli=x.col(i);
							stable_sort(coli.begin(),coli.end(),greater<double>());
							f.col(i)=coli;
					    }
					}else{
						colvec coli(n);
						for(int i=0;i<p;++i){
							coli=x.col(i);
							stable_sort(coli.begin(),coli.end(),greater<double>());
							f.col(i)=coli;
					    }
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<p;++i){
							colvec coli=x.col(i);
							sort(coli.begin(),coli.end(),greater<double>());
							f.col(i)=coli;
					    }
					}else{
						colvec coli(n);
						for(int i=0;i<p;++i){
							coli=x.col(i);
							sort(coli.begin(),coli.end(),greater<double>());
							f.col(i)=coli;
					    }
					}
				}
			}
			else{
				if(stable){
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<p;++i){
							colvec coli=x.col(i);
							stable_sort(coli.begin(),coli.end());
							f.col(i)=coli;
					    }
					}else{
						colvec coli(n);
						for(int i=0;i<p;++i){
							coli=x.col(i);
							stable_sort(coli.begin(),coli.end());
							f.col(i)=coli;
					    }
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<p;++i){
							colvec coli=x.col(i);
							sort(coli.begin(),coli.end());
							f.col(i)=coli;
					    }
					}else{
						colvec coli(n);
						for(int i=0;i<p;++i){
							coli=x.col(i);
							sort(coli.begin(),coli.end());
							f.col(i)=coli;
					    }
					}
				}
			}  
			return f;
		}

		inline NumericMatrix rowSort(NumericMatrix x,const bool descend=false,const bool stable=false,const bool parallel=false){
			const int n=x.nrow(),p=x.ncol();
			NumericMatrix f(n,p);
			mat xx(x.begin(),n,p,false),ff(f.begin(),n,p,false);
			if(descend){
				if(stable){
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<n;++i){
							rowvec rowi=xx.row(i);
							stable_sort(rowi.begin(),rowi.end(),greater<double>());
							ff.row(i)=rowi;
					    }
					}else{
						rowvec rowi(n);
						for(int i=0;i<n;++i){
							rowi=xx.row(i);
							stable_sort(rowi.begin(),rowi.end(),greater<double>());
							ff.row(i)=rowi;
					    }
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<n;++i){
							rowvec rowi=xx.row(i);
							sort(rowi.begin(),rowi.end(),greater<double>());
							ff.row(i)=rowi;
					    }
					}else{
						rowvec rowi(n);
						for(int i=0;i<n;++i){
							rowi=xx.row(i);
							sort(rowi.begin(),rowi.end(),greater<double>());
							ff.row(i)=rowi;
					    }
					}
				}
			}
			else{
				if(stable){
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<n;++i){
							rowvec rowi=xx.row(i);
							stable_sort(rowi.begin(),rowi.end());
							ff.row(i)=rowi;
					    }
					}else{
						rowvec rowi(n);
						for(int i=0;i<n;++i){
							rowi=xx.row(i);
							stable_sort(rowi.begin(),rowi.end());
							ff.row(i)=rowi;
					    }
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<n;++i){
							rowvec rowi=xx.row(i);
							sort(rowi.begin(),rowi.end());
							ff.row(i)=rowi;
					    }
					}else{
						rowvec rowi(n);
						for(int i=0;i<n;++i){
							rowi=xx.row(i);
							sort(rowi.begin(),rowi.end());
							ff.row(i)=rowi;
					    }
					}
				}
			}  
			return f;
		}

		inline mat rowSort(mat x,const bool descend=false,const bool stable=false,const bool parallel=false){
			const int n=x.n_rows,p=x.n_cols;
			mat f(n,p);
			if(descend){
				if(stable){
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<n;++i){
							rowvec rowi=x.row(i);
							stable_sort(rowi.begin(),rowi.end(),greater<double>());
							f.row(i)=rowi;
					    }
					}else{
						rowvec rowi(n);
						for(int i=0;i<n;++i){
							rowi=x.row(i);
							stable_sort(rowi.begin(),rowi.end(),greater<double>());
							f.row(i)=rowi;
					    }
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<n;++i){
							rowvec rowi=x.row(i);
							sort(rowi.begin(),rowi.end(),greater<double>());
							f.row(i)=rowi;
					    }
					}else{
						rowvec rowi(n);
						for(int i=0;i<n;++i){
							rowi=x.row(i);
							sort(rowi.begin(),rowi.end(),greater<double>());
							f.row(i)=rowi;
					    }
					}
				}
			}
			else{
				if(stable){
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<n;++i){
							rowvec rowi=x.row(i);
							stable_sort(rowi.begin(),rowi.end());
							f.row(i)=rowi;
					    }
					}else{
						rowvec rowi(n);
						for(int i=0;i<n;++i){
							rowi=x.row(i);
							stable_sort(rowi.begin(),rowi.end());
							f.row(i)=rowi;
					    }
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
					    #pragma omp parallel for
					    #endif
						for(int i=0;i<n;++i){
							rowvec rowi=x.row(i);
							sort(rowi.begin(),rowi.end());
							f.row(i)=rowi;
					    }
					}else{
						rowvec rowi(n);
						for(int i=0;i<n;++i){
							rowi=x.row(i);
							sort(rowi.begin(),rowi.end());
							f.row(i)=rowi;
					    }
					}
				}
			}  
			return f;
		}

		inline bool is_symmetric(NumericMatrix x){
			int ncl=x.ncol(),i,j;
			for(i=1;i<ncl;++i){
				for(j=0;j<i;++j){
					if(x(j,i)!=x(i,j)){
						return false;
					}
				}
			}
			return true;
		}

		inline bool is_symmetric(mat x){
			int ncl=x.n_cols,i,j;
			for(i=1;i<ncl;++i){
				for(j=0;j<i;++j){
					if(x(j,i)!=x(i,j)){
						return false;
					}
				}
			}
			return true;
		}

		inline NumericVector colMedian(NumericMatrix x,const bool na_rm=false,const bool parallel=false){
		    return na_rm ? col_meds_na_rm(x,parallel) : col_meds_simple(x,parallel);
		}

		inline rowvec colMedian(mat x,const bool na_rm=false,const bool parallel=false){
		    return na_rm ? col_meds_na_rm(x,parallel) : col_meds_simple(x,parallel);
		}

		inline NumericVector rowMedian(NumericMatrix x,const bool na_rm=false,const bool parallel=false){
		    return na_rm ? row_meds_na_rm(x,parallel) : row_meds_simple(x,parallel);
		}

		inline colvec rowMedian(mat x,const bool na_rm=false,const bool parallel=false){
		    return na_rm ? row_meds_na_rm(x,parallel) : row_meds_simple(x,parallel);
		}

		inline NumericVector colVars(NumericMatrix x,const bool std=false,const bool na_rm=false,const bool parallel=false){
		    mat xx(x.begin(),x.nrow(),x.ncol(),false);
		    NumericVector f(xx.n_cols);
		    rowvec ff(f.begin(),f.size(),false);
		    if(parallel){
				#pragma omp parallel for
		        for(unsigned int i=0;i<xx.n_cols;++i){
		            ff[i]=Rfast::vector::var<colvec>(xx.col(i),std,na_rm);
		        }
		    }else{
		        for(unsigned int i=0;i<xx.n_cols;++i){
		            ff[i]=Rfast::vector::var<colvec>(xx.col(i),std,na_rm);
		        }
		    }
		    return f;
		}

		inline NumericVector rowVars(NumericMatrix x,const bool std=false,const bool na_rm=false,const bool parallel=false){
		    mat xx(x.begin(),x.nrow(),x.ncol(),false);
		    NumericVector f(xx.n_rows);
		    colvec ff(f.begin(),f.size(),false);
		    if(parallel){
		    	#ifdef _OPENMP
				#pragma omp parallel for
				#endif
		        for(unsigned int i=0;i<xx.n_rows;++i){
		            ff[i]=Rfast::vector::var<rowvec>(xx.row(i),std,na_rm);
		        }
		    }else{
		        for(unsigned int i=0;i<xx.n_rows;++i){
		            ff[i]=Rfast::vector::var<rowvec>(xx.row(i),std,na_rm);
		        }
		    }
		    return f;
		}

		inline rowvec colVars(mat x,const bool std=false,const bool na_rm=false,const bool parallel=false){
		    rowvec f(x.n_cols);
		    if(parallel){
		    	#ifdef _OPENMP
				#pragma omp parallel for
				#endif
		        for(unsigned int i=0;i<x.n_cols;++i){
		            f[i]=Rfast::vector::var<colvec>(x.col(i),std,na_rm);
		        }
		    }else{
		        for(unsigned int i=0;i<x.n_cols;++i){
		            f[i]=Rfast::vector::var<colvec>(x.col(i),std,na_rm);
		        }
		    }
		    return f;
		}

		inline colvec rowVars(mat x,const bool std=false,const bool na_rm=false,const bool parallel=false){
		    colvec f(x.n_rows);
		    if(parallel){
				#pragma omp parallel for
		        for(unsigned int i=0;i<x.n_rows;++i){
		            f[i]=Rfast::vector::var<rowvec>(x.row(i),std,na_rm);
		        }
		    }else{
		        for(unsigned int i=0;i<x.n_rows;++i){
		            f[i]=Rfast::vector::var<rowvec>(x.row(i),std,na_rm);
		        }
		    }
		    return f;
		}

		inline NumericVector colMads(NumericMatrix x,const string method="median",const bool na_rm=false,const bool parallel=false){
		    mat xx(x.begin(),x.nrow(),x.ncol(),false);
		    NumericVector f(xx.n_cols);
		    rowvec ff(f.begin(),f.size(),false);
		    if(parallel){
		    #pragma omp parallel for
		        for(unsigned int i=0;i<xx.n_cols;++i){
		            ff[i]=Rfast::vector::mad<colvec>(xx.col(i),method,na_rm);
		        }
		    }else{
		        for(unsigned int i=0;i<xx.n_cols;++i){
		            ff[i]=Rfast::vector::mad<colvec>(xx.col(i),method,na_rm);
		        }
		    }
		    return f;
		}


		inline NumericVector rowMads(NumericMatrix x,const string method="median",const bool na_rm=false,const bool parallel=false){
		    mat xx(x.begin(),x.nrow(),x.ncol(),false);
		    NumericVector f(xx.n_rows);
		    colvec ff(f.begin(),f.size(),false);
		    if(parallel){
		        #pragma omp parallel for
		        for(unsigned int i=0;i<xx.n_rows;++i){
		            ff[i]=Rfast::vector::mad<rowvec>(xx.row(i),method,na_rm);
		        }
		    }else{
		        for(unsigned int i=0;i<xx.n_rows;++i){
		            ff[i]=Rfast::vector::mad<rowvec>(xx.row(i),method,na_rm);
		        }
		    }
		    return f;
		}

		inline rowvec colMads(mat x,const string method="median",const bool na_rm=false,const bool parallel=false){
		    rowvec f(x.n_cols);
		    if(parallel){
		    #pragma omp parallel for
		        for(unsigned int i=0;i<x.n_cols;++i){
		            f[i]=Rfast::vector::mad<colvec>(x.col(i),method,na_rm);
		        }
		    }else{
		        for(unsigned int i=0;i<x.n_cols;++i){
		            f[i]=Rfast::vector::mad<colvec>(x.col(i),method,na_rm);
		        }
		    }
		    return f;
		}

		inline colvec rowMads(mat x,const string method="median",const bool na_rm=false,const bool parallel=false){
		    colvec f(x.n_rows);
		    if(parallel){
		        #pragma omp parallel for
		        for(unsigned int i=0;i<x.n_rows;++i){
		            f[i]=Rfast::vector::mad<rowvec>(x.row(i),method,na_rm);
		        }
		    }else{
		        for(unsigned int i=0;i<x.n_rows;++i){
		            f[i]=Rfast::vector::mad<rowvec>(x.row(i),method,na_rm);
		        }
		    }
		    return f;
		}

		template<class Engine=std::default_random_engine>
		inline NumericMatrix colShuffle(NumericMatrix x,Engine engine=Engine()){
		    const int n=x.ncol();
		    seed_seq seq{get_current_nanoseconds()};
		    std::vector<long long unsigned  int> seeds(n);
		    seq.generate(seeds.begin(),seeds.end());
		    NumericMatrix y(x.nrow(),n);
		    for(int i=0;i<n;++i){
		        engine.seed(seeds[i]);
		        y.column(i)=Rfast::vector::shuffle<NumericVector>(x.column(i),engine);
		    }
		    return y;
		}

		template<class Engine=std::default_random_engine>
		NumericMatrix rowShuffle(NumericMatrix x,Engine engine=Engine()){
		    const int n=x.ncol();
		    seed_seq seq{get_current_nanoseconds()};
		    std::vector<long long unsigned  int> seeds(n);
		    seq.generate(seeds.begin(),seeds.end());
		    NumericMatrix y(n,x.ncol());
		    for(int i=0;i<n;++i){
		        engine.seed(seeds[i]);
		        y.row(i)=Rfast::vector::shuffle<NumericVector>(x.row(i),engine);
		    }
		    return y;
		}

		template<class Engine=std::default_random_engine>
		mat colShuffle(mat x,Engine engine=Engine()){
		    const int n=x.n_cols;
		    seed_seq seq{get_current_nanoseconds()};
		    std::vector<long long unsigned int> seeds(n);
		    seq.generate(seeds.begin(),seeds.end());
		    mat y(x.n_rows,n);
		    for(int i=0;i<n;++i){
		        engine.seed(seeds[i]);
		        y.col(i)=Rfast::vector::shuffle<colvec>(x.col(i),engine);
		    }
		    return y;
		}

		template<class Engine=std::default_random_engine>
		mat rowShuffle(mat x,Engine engine=Engine()){
		    const int n=x.n_rows;
		    seed_seq seq{get_current_nanoseconds()};
		    std::vector<long long unsigned int> seeds(n);
		    seq.generate(seeds.begin(),seeds.end());
		    mat y(n,x.n_cols);
		    for(int i=0;i<n;++i){
				engine.seed(seeds[i]);
		        y.row(i)=Rfast::vector::shuffle<rowvec>(x.row(i),engine);
		    }
		    return y;
		}
	}
}

#endif
