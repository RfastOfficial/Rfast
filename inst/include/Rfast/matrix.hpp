
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
#include "types.hpp"
#include "helpers.hpp"
#include "vector.hpp"

namespace Rfast {
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
			int p=0,n=0;

			if(tx){
				n=x.n_cols;
				p=x.n_cols;
			}else{
				n=x.n_rows;
				p= ty ? x.n_rows : x.n_cols;
			}
			mat C(n,p);
			colvec yi;

			if(!tx and !ty){ // matrix multiplication
				mat xx=Rfast::transpose(x);
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
				mat yy=Rfast::transpose(y);
				mat xx=Rfast::transpose(x);
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

		inline 
		DataFrame colSort(DataFrame x,const bool descend,const bool stable, const bool parallel){
			List f(x.size());
			if(descend){
				if(stable){
					if(parallel){
						#ifdef _OPENMP
						#pragma omp parallel for
						#endif
						for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
							switch(Type::type(s->get())){
								case Type::Types::REAL:
								setResultParallelSection<colvec,NumericVector,std::stable_sort>(f,s,mgreater<bool,double,double>);
								break;
								case Type::Types::INT:
								setResultParallelSection<icolvec,IntegerVector,std::stable_sort>(f,s,mgreater<bool,int,int>);
								break;
								case Type::Types::CHAR:
								setResultParallelSection<vector<string>,CharacterVector,std::stable_sort>(f,s,mgreater<bool,string,string>);
								break;
								case Type::Types::FACTOR:
								setResultParallelSection<vector<string>,CharacterVector,std::stable_sort>(f,s,mgreater<bool,string,string>);
								break;
								default:
									break;
							}
						}
					}else{
						int i=0;
						for(auto c : x){
							switch(Type::type(c)){
								case Type::Types::REAL:
								setResult<NumericVector,std::stable_sort>(f,i++,c,mgreater<bool,double,double>);
								break;
								case Type::Types::INT:
								setResult<IntegerVector,std::stable_sort>(f,i++,c,mgreater<bool,int,int>);
								break;
								case Type::Types::CHAR:
								setResult<CharacterVector,std::stable_sort>(f,i++,c,mgreater<bool,CharacterVector::value_type,CharacterVector::value_type>);
								break;
								case Type::Types::FACTOR:
								setResult<CharacterVector,std::stable_sort>(f,i++,c,mgreater<bool,CharacterVector::value_type,CharacterVector::value_type>);
								break;
								default:
									break;
							}
						}
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
						#pragma omp parallel for
						#endif
						for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
							switch(Type::type(s->get())){
								case Type::Types::REAL:
								setResultParallelSection<colvec,NumericVector,std::sort>(f,s,mgreater<bool,double,double>);
								break;
								case Type::Types::INT:
								setResultParallelSection<icolvec,IntegerVector,std::sort>(f,s,mgreater<bool,int,int>);
								break;
								case Type::Types::CHAR:
								setResultParallelSection<vector<string>,CharacterVector,std::sort>(f,s,mgreater<bool,string,string>);
								break;
								case Type::Types::FACTOR:
								setResultParallelSection<vector<string>,CharacterVector,std::sort>(f,s,mgreater<bool,string,string>);
								break;
								default:
									break;
							}
						}
					}else{
						int i=0;
						for(auto c : x){
							switch(Type::type(c)){
								case Type::Types::REAL:
								setResult<NumericVector,std::sort>(f,i++,c,mgreater<bool,double,double>);
								break;
								case Type::Types::INT:
								setResult<IntegerVector,std::sort>(f,i++,c,mgreater<bool,int,int>);
								break;
								case Type::Types::CHAR:
								setResult<CharacterVector,std::sort>(f,i++,c,mgreater<bool,CharacterVector::value_type,CharacterVector::value_type>);
								break;
								case Type::Types::FACTOR:
								setResult<CharacterVector,std::sort>(f,i++,c,mgreater<bool,CharacterVector::value_type,CharacterVector::value_type>);
								break;
								default:
									break;
							}
						}
					}
				}
			}else{
				if(stable){
					if(parallel){
						#ifdef _OPENMP
						#pragma omp parallel for
						#endif
						for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
							switch(Type::type(s->get())){
								case Type::Types::REAL:
								setResultParallelSection<colvec,NumericVector,std::stable_sort>(f,s);
								break;
								case Type::Types::INT:
								setResultParallelSection<icolvec,IntegerVector,std::stable_sort>(f,s);
								break;
								case Type::Types::CHAR:
								setResultParallelSection<vector<string>,CharacterVector,std::stable_sort>(f,s);
								break;
								case Type::Types::FACTOR:
								setResultParallelSection<vector<string>,CharacterVector,std::stable_sort>(f,s);
								break;
								default:
									break;
							}
						}
					}else{
						int i=0;
						for(auto c : x){
							switch(Type::type(c)){
								case Type::Types::REAL:
								setResult<NumericVector,std::stable_sort>(f,i++,c);
								break;
								case Type::Types::INT:
								setResult<IntegerVector,std::stable_sort>(f,i++,c);
								break;
								case Type::Types::CHAR:
								setResult<CharacterVector,std::stable_sort>(f,i++,c);
								break;
								case Type::Types::FACTOR:
								setResult<CharacterVector,std::stable_sort>(f,i++,c);
								break;
								default:
									break;
							}
						}
					}
				}else{
					if(parallel){
						#ifdef _OPENMP
						#pragma omp parallel for
						#endif
						for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
							switch(Type::type(s->get())){
								case Type::Types::REAL:
								setResultParallelSection<colvec,NumericVector,std::sort>(f,s);
								break;
								case Type::Types::INT:
								setResultParallelSection<icolvec,IntegerVector,std::sort>(f,s);
								break;
								case Type::Types::CHAR:
								setResultParallelSection<vector<string>,CharacterVector,std::sort>(f,s);
								break;
								case Type::Types::FACTOR:
								setResultParallelSection<vector<string>,CharacterVector,std::sort>(f,s);
								break;
								default:
									break;
							}
						}
					}else{
						int i=0;
						for(auto c : x){
							switch(Type::type(c)){
								case Type::Types::REAL:
								setResult<NumericVector,std::sort>(f,i++,c);
								break;
								case Type::Types::INT:
								setResult<IntegerVector,std::sort>(f,i++,c);
								break;
								case Type::Types::CHAR:
								setResult<CharacterVector,std::sort>(f,i++,c);
								break;
								case Type::Types::FACTOR:
								setResult<CharacterVector,std::sort>(f,i++,c);
								break;
								default:
									break;
							}
						}
					}
				}
			}
			f.names() = x.names();
			return f;
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

		inline NumericVector colMedian(DataFrame& x,const bool na_rm=false,const bool parallel=false){
			NumericVector f(x.size());
			if(na_rm){
				if(parallel){
					colvec ff(f.begin(),f.size(),false);
					#pragma omp parallel for
					for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
						colvec y;
						int i;
						#pragma omp critical
						{
							NumericVector yy;
							yy=*s;
							y = colvec(yy.begin(),yy.size());
							i = s-x.begin();
						}
						ff[i]=med_helper<colvec>(y.begin(),y.begin()+(int)(std::remove_if(y.begin(),y.end(),R_IsNA)-y.begin()));
					}
				}else{
					int i=0;
					NumericVector y(x.nrows());
					for(auto c : x){
						y=c;
						f[i++]=med_helper<NumericVector>(y.begin(),y.begin()+(int)(std::remove_if(y.begin(),y.end(),R_IsNA)-y.begin()));
					}
				}
			}else{
				const int step=x.nrow(),middle=step/2-1;
				if(step%2==0){
					if(parallel){
						colvec ff(f.begin(),f.size(),false);
						#pragma omp parallel for
						for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
							colvec y;
							int i;
							#pragma omp critical
							{
								NumericVector yy;
								yy=*s;
								y = colvec(yy.begin(),yy.size());
								i = s-x.begin();
							}
							nth_element(y.begin(),y.begin()+middle,y.end());
							ff[i]=(y[middle]+*(min_element(y.begin()+middle+1,y.end())))/2.0;
						}
					}else{
						int i=0;
						NumericVector y(x.nrows());
						for(auto c : x){
							y=c;
							nth_element(y.begin(),y.begin()+middle,y.end());
							f[i++]=(y[middle]+*(min_element(y.begin()+middle+1,y.end())))/2.0;
						}
					}
				}else{
					if(parallel){
						colvec ff(f.begin(),f.size(),false);
						#pragma omp parallel for
						for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
							colvec y;
							int i;
							#pragma omp critical
							{
								NumericVector yy;
								yy=*s;
								y = colvec(yy.begin(),yy.size());
								i = s-x.begin();
							}
							nth_element(y.begin(),y.begin()+middle+1,y.end());
							ff[i]=y[middle+1];
						}
					}else{
						int i=0;
						NumericVector y(x.nrows());
						for(auto c : x){
							y=c;
							nth_element(y.begin(),y.begin()+middle+1,y.end());
							f[i++]=NumericVector(y.begin(),y.end())[middle+1];
						}
					}
				}
			}
			f.names() = x.names();
			return f;
		}

		inline NumericVector colMedian(NumericMatrix& x,const bool na_rm=false,const bool parallel=false){
			const int p=x.ncol();
			NumericVector F(p);
			if(na_rm){
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
			}else{
				const int step=x.nrow(),middle=step/2-1;
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
			}
			return F;
		}

		inline rowvec colMedian(mat& x,const bool na_rm=false,const bool parallel=false){
			const int p=x.n_cols;
			rowvec F(p);
			if(na_rm){
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
			}else{
				const int step=x.n_rows,middle=step/2-1;
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
			}
			return F;
		}

		inline NumericVector rowMedian(NumericMatrix x,const bool na_rm=false,const bool parallel=false){
			const int p=x.nrow();
			NumericVector F(p);
			if(na_rm){
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
			}else{
				const int sz=x.ncol(),middle=sz/2-1;
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
			}
			return F;
		}

		inline colvec rowMedian(mat& x,const bool na_rm=false,const bool parallel=false){
			const int p=x.n_rows;
			colvec F(p);
			if(na_rm){
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
			}else{
				const int sz=x.n_cols,middle=sz/2-1;
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
			}
			return F;
		}

		inline NumericVector colVars(NumericMatrix x,const bool std=false,const bool na_rm=false,const bool parallel=false){
			mat xx(x.begin(),x.nrow(),x.ncol(),false);
			NumericVector f(xx.n_cols);
			rowvec ff(f.begin(),f.size(),false);
			if(parallel){
				#pragma omp parallel for
				for(unsigned int i=0;i<xx.n_cols;++i){
					ff[i]=Rfast::var<colvec>(xx.col(i),std,na_rm);
				}
			}else{
				for(unsigned int i=0;i<xx.n_cols;++i){
					ff[i]=Rfast::var<colvec>(xx.col(i),std,na_rm);
				}
			}
			return f;
		}

		inline NumericVector colVars(DataFrame x,const bool std=false,const bool na_rm=false,const bool parallel=false){
			NumericVector f(x.size());
			if(parallel){
				colvec ff(f.begin(),f.size(),false);
				#pragma omp parallel for
				for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
					colvec y;
					int i;
					#pragma omp critical
					{
						NumericVector yy;
						yy=*s;
						y = colvec(yy.begin(),yy.size(),false);
						i = s-x.begin();
					}
					ff[i]=Rfast::var<colvec>(y,std,na_rm);
				}
			}else{
				int i=0;
				NumericVector y(x.nrows());
				for(auto c : x){
					y=c;
					f[i++]=Rfast::var<colvec>(y,std,na_rm);
				}
			}
			f.names() = x.names();
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
					ff[i]=Rfast::var<rowvec>(xx.row(i),std,na_rm);
				}
			}else{
				for(unsigned int i=0;i<xx.n_rows;++i){
					ff[i]=Rfast::var<rowvec>(xx.row(i),std,na_rm);
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
					f[i]=Rfast::var<colvec>(x.col(i),std,na_rm);
				}
			}else{
				for(unsigned int i=0;i<x.n_cols;++i){
					f[i]=Rfast::var<colvec>(x.col(i),std,na_rm);
				}
			}
			return f;
		}

		inline colvec rowVars(mat x,const bool std=false,const bool na_rm=false,const bool parallel=false){
			colvec f(x.n_rows);
			if(parallel){
				#pragma omp parallel for
				for(unsigned int i=0;i<x.n_rows;++i){
					f[i]=Rfast::var<rowvec>(x.row(i),std,na_rm);
				}
			}else{
				for(unsigned int i=0;i<x.n_rows;++i){
					f[i]=Rfast::var<rowvec>(x.row(i),std,na_rm);
				}
			}
			return f;
		}

		inline NumericVector colMads(DataFrame x,const string method="median",const bool na_rm=false,const bool parallel=false){
			NumericVector f(x.size());
			if(parallel){
				colvec ff(f.begin(),f.size(),false);
				#pragma omp parallel for
				for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
					colvec y;
					int i;
					#pragma omp critical
					{
						NumericVector yy;
						yy=*s;
						y = colvec(yy.begin(),yy.size(),false);
						i = s-x.begin();
					}
					ff[i]=Rfast::mad<colvec>(y,method,na_rm);
				}
			}else{
				int i=0;
				NumericVector y(x.nrows());
				for(auto c : x){
					y=c;
					f[i++]=Rfast::mad<colvec>(y,method,na_rm);
				}
			}
			f.names() = x.names();
			return f;
		}

		inline NumericVector colMads(NumericMatrix x,const string method="median",const bool na_rm=false,const bool parallel=false){
			mat xx(x.begin(),x.nrow(),x.ncol(),false);
			NumericVector f(xx.n_cols);
			rowvec ff(f.begin(),f.size(),false);
			if(parallel){
		    #pragma omp parallel for
				for(unsigned int i=0;i<xx.n_cols;++i){
					ff[i]=Rfast::mad<colvec>(xx.col(i),method,na_rm);
				}
			}else{
				for(unsigned int i=0;i<xx.n_cols;++i){
					ff[i]=Rfast::mad<colvec>(xx.col(i),method,na_rm);
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
					ff[i]=Rfast::mad<rowvec>(xx.row(i),method,na_rm);
				}
			}else{
				for(unsigned int i=0;i<xx.n_rows;++i){
					ff[i]=Rfast::mad<rowvec>(xx.row(i),method,na_rm);
				}
			}
			return f;
		}

		inline rowvec colMads(mat x,const string method="median",const bool na_rm=false,const bool parallel=false){
			rowvec f(x.n_cols);
			if(parallel){
		    #pragma omp parallel for
				for(unsigned int i=0;i<x.n_cols;++i){
					f[i]=Rfast::mad<colvec>(x.col(i),method,na_rm);
				}
			}else{
				for(unsigned int i=0;i<x.n_cols;++i){
					f[i]=Rfast::mad<colvec>(x.col(i),method,na_rm);
				}
			}
			return f;
		}

		inline colvec rowMads(mat x,const string method="median",const bool na_rm=false,const bool parallel=false){
			colvec f(x.n_rows);
			if(parallel){
		        #pragma omp parallel for
				for(unsigned int i=0;i<x.n_rows;++i){
					f[i]=Rfast::mad<rowvec>(x.row(i),method,na_rm);
				}
			}else{
				for(unsigned int i=0;i<x.n_rows;++i){
					f[i]=Rfast::mad<rowvec>(x.row(i),method,na_rm);
				}
			}
			return f;
		}

	template<class Engine=std::default_random_engine>
		inline DataFrame colShuffle(DataFrame x,Engine engine=Engine()){
			const int n = x.size();
			seed_seq seq{get_current_nanoseconds()};
			std::vector<long long unsigned  int> seeds(n);
			seq.generate(seeds.begin(),seeds.end());
			List f(n);
			// if(parallel){
			// 	colvec ff(f.begin(),f.size(),false);
			// 	#pragma omp parallel for
			// 	for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
			// 		colvec y;
			// 		int i;
			// 		#pragma omp critical
			// 		{
			// 			engine.seed(seeds[i]);
			// 			NumericVector yy;
			// 			yy=*s;
			// 			y = colvec(yy.begin(),yy.size(),false);
			// 			i = s-x.begin();
			// 		}
			// 		ff[i]=Rfast::shuffle<NumericVector>(y,engine);
			// 	}
			// }else{
			int i=0;
			NumericVector y(x.nrows());
			for(auto c : x){
				y=c;
				engine.seed(seeds[i]);
				f[i++]=Rfast::shuffle<NumericVector>(y,engine);
			}
			// }
			f.names() = x.names();
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
				y.column(i)=Rfast::shuffle<NumericVector>(x.column(i),engine);
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
				y.row(i)=Rfast::shuffle<NumericVector>(x.row(i),engine);
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
				y.col(i)=Rfast::shuffle<colvec>(x.col(i),engine);
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
				y.row(i)=Rfast::shuffle<rowvec>(x.row(i),engine);
			}
			return y;
		}
	}

#endif
