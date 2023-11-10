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
            DATAFRAME,
            LOGICAL
        };

        struct R {
            inline static const int Real = REALSXP;
            inline static const int Int = INTSXP;
            inline static const int Char = CHARSXP;
            inline static const int Lgl = LGLSXP;
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
                case LGLSXP: return Types::LOGICAL;
                default:
                stop("Error: unsupported type.\n");
            }
        }
    };
    namespace internal { // the initialization of static variables are defined in file types.cpp
	    template<class T>
	    struct NA_helper : std::false_type {
	    	using type = T;
	    	inline static type val;
	    };
	    template<>
	    struct NA_helper<int> : std::true_type {
	    	using type = int;
	    	inline static type val= NA_INTEGER;
	    };
	    template<>
	    struct NA_helper<double> : std::true_type {
	    	using type = double;
	    	inline static type val= NA_LOGICAL;
	    };
	    template<>
	    struct NA_helper<bool> : std::true_type {
	    	using type = int;
	    	inline static type val= NA_REAL;
	    };
	    template<>
	    struct NA_helper<string> : std::true_type {
	    	using type = SEXP;
	    	inline static type val = NA_STRING;
	    };
	}

    template<class T>
    struct NA {
        static typename internal::NA_helper<T>::type value(){
        	static_assert(internal::NA_helper<T>{}, "Unsupported type for NA.");
    		return internal::NA_helper<T>::val;
        }
    };

    // struct Na {
    //     inline operator typename internal::NA_helper<double>::type(){
    //         return internal::NA_helper<double>::val;
    //     }
    //     inline operator typename internal::NA_helper<string>::type(){
    //         return internal::NA_helper<string>::val;
    //     }
    //     inline operator typename internal::NA_helper<bool>::type(){
    //         return internal::NA_helper<bool>::val;
    //     }
    //     inline operator typename internal::NA_helper<double>::type(){
    //         return internal::NA_helper<double>::val;
    //     }
    // };
}

#endif