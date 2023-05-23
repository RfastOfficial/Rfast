#ifndef TYPES_H
#define TYPES_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <type_traits>
using namespace Rcpp;
using namespace arma;

namespace Rfast {
    struct Type
    {
        enum class Types {
            REAL,
            INT,
            CHAR,
            FACTOR,
            LIST,
            DATAFRAME
        };

        template<class T>
        static Types type(T t){
            if(Rf_isFactor(t)) return Types::FACTOR;
            if(Rf_isNewList(t)) return Types::DATAFRAME;
            switch(TYPEOF(t)){
                case REALSXP: return Types::REAL;
                case INTSXP: return Types::INT;
                case CHARSXP: return Types::CHAR;
                case LISTSXP: return Types::LIST;
                default:
                stop("Error: unsupported type.\n");
            }
        }
    };
    namespace internal { // the initialization of static variables are defined in file types.cpp
	    template<class T>
	    struct NA_helper : std::false_type {

	    	using type = T;
	    	static type val;
	    };
	    template<>
	    struct NA_helper<int> : std::true_type {
	    	using type = int;
	    	static type val;
	    };
	    template<>
	    struct NA_helper<double> : std::true_type {
	    	using type = double;
	    	static type val;
	    };
	    template<>
	    struct NA_helper<bool> : std::true_type {
	    	using type = int;
	    	static type val;
	    };
	    template<>
	    struct NA_helper<string> : std::true_type {
	    	using type = SEXP;
	    	static type val;
	    };
	}

    template<class T>
    struct NA {
        static typename internal::NA_helper<T>::type value(){
        	static_assert(internal::NA_helper<T>{}, "Unsupported type for NA.");
    		return internal::NA_helper<T>::val;
        }
    };

    
    inline Rfast::internal::NA_helper<string>::type Rfast::internal::NA_helper<string>::val = NA_STRING;
    inline Rfast::internal::NA_helper<int>::type    Rfast::internal::NA_helper<int>::val    = NA_INTEGER;
    inline Rfast::internal::NA_helper<bool>::type   Rfast::internal::NA_helper<bool>::val   = NA_LOGICAL;
    inline Rfast::internal::NA_helper<double>::type Rfast::internal::NA_helper<double>::val = NA_REAL;
}

#endif