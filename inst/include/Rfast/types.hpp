#ifndef TYPES_H
#define TYPES_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <type_traits>
#include <string>
using std::string;
using namespace Rcpp;
using namespace arma;

namespace Rfast {

    struct R {
        inline static const int Real = REALSXP;
        inline static const int Int = INTSXP;
        inline static const int Char = CHARSXP;
        inline static const int String = STRSXP;
        inline static const int List = LISTSXP;
        inline static const int Lgl = LGLSXP;
        inline static const int Complex = CPLXSXP;
        inline static SEXP Null = R_NilValue;
    };

    namespace Type
    {
        enum class Types {
            REAL,
            INT,
            CHAR,
            STRING,
            FACTOR,
            LIST,
            DATAFRAME,
            LOGICAL,
            COMPLEX
        };

        template<class T, class U>
        static Types type([[maybe_unused]] U t){
            if constexpr(std::is_same<T, int>::value){
                return Types::INT;
            }else if constexpr(std::is_same<T, double>::value){
                return Types::REAL;
            }else if constexpr(std::is_same<T, string>::value){
                return Types::STRING;
            }else if constexpr(std::is_same<T, bool>::value){
                return Types::LOGICAL;
            }else if constexpr(std::is_same<T, char>::value){
                return Types::CHAR;
            }else{
                stop("Error: unsupported type.\n");
            }
        }
        
        template<> [[maybe_unused]]
        Types type<SEXP>(SEXP t){
            if(Rf_isFactor(t)) return Types::FACTOR;
            if(Rf_isNewList(t)) return Types::DATAFRAME;
            switch(TYPEOF(t)){
                case R::Real: return Types::REAL;
                case R::Int: return Types::INT;
                case R::Char: return Types::CHAR;
                case R::String: return Types::STRING;
                case R::List: return Types::LIST;
                case R::Lgl: return Types::LOGICAL;
                case R::Complex: return Types::COMPLEX;
                default:
                stop("Error: unsupported type.\n");
            }
        }
    };

    template<class T> 
    typename std::conditional<std::is_same<T, SEXP>::value, SEXP, T*>::type asPtr(SEXP x){
        if constexpr(std::is_same<T,SEXP>::value){
            return x;
        }else if constexpr(std::is_same<T,double>::value){
            return REAL(x);
        }else if constexpr(std::is_same<T,int>::value){
            return INTEGER(x);
        }else if constexpr(std::is_same<T,Rcomplex>::value){
            return COMPLEX(x);
        }
        return nullptr;
    }

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