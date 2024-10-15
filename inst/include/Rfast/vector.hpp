
#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "templates.h"

namespace Rfast {

	using namespace Rcpp;
	using namespace arma;
	using std::remove_if;
	using std::string;

		//type T can any iterable class
	template<typename T>
	inline double median(T x){
		return med_helper<T>(x.begin(),x.end());
	}

		//type T can be any vector
	template<class T>
	inline double var(T x,const bool std=false,const bool na_rm=false){
		int n = 0;
		double v=0.0,sum1=0,sum2=0;
		if(na_rm){
			for(auto v : x){
				if(!R_IsNA(v)){
					sum1+=v*v;
					sum2+=v;
					++n;
				}
			}
		}else{
			n=x.size();
			double *xx=&x[0],*end=xx+n;
			for(;end-xx;++xx){
				v=*xx;
				sum1+=v*v;
				sum2+=v;
			}
		}
		v=(sum1-sum2*sum2/n)/(n-1);
		return std ? sqrt(v) : v;
	}

		// type T can any iterable class
		// argument method can be one of "median","mean"
		// argument na_rm for remove NAs from the vector using R's "R_isNA" function
	template<typename T>
	inline double mad(T x,string method="median",const bool na_rm=false){
		const int newsize = na_rm ? remove_if(x.begin(),x.end(),R_IsNA)-x.begin() : x.size();
		double res=0;
		if(newsize>1){
			T xx(x.begin(),newsize,false);
			if(method=="median"){
				const double center = 1.4826;
				const double md=med_helper<T>(xx.begin(),xx.end());
				T y=abs(xx-md);
				res=med_helper<T>(y.begin(),y.end())*center;
			}else if(method=="mean"){
				res=mean(abs(xx-mean(xx)));
			}else{
				stop("Wrong method. Choose \"median\" or \"mean\"");
			}
		}else{
			res=NA<double>::value();
		}
		return res;
	}
		// type T can any iterable class
		// argument method can be one of "median","mean"
		// argument na_rm for remove NAs from the vector using R's "R_isNA" function
	template<>
	inline double mad<NumericVector>(NumericVector xx ,const string method,const bool na_rm){
		colvec x(xx.begin(),newsize,false);
		return mad<colvec>(x,method,na_rm);
	}
		// type T can any iterable class
		// type Engine can by anyone that overloads operator()
	template<class T,class Engine=std::default_random_engine>
	inline T shuffle(T x,Engine engine=Engine(get_current_nanoseconds())){
		std::shuffle(x.begin(),x.end(),engine);
		return x;
	}
}
#endif