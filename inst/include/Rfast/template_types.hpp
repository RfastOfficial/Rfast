
#pragma once

#include <RcppArmadillo.h>


template<typename f,typename s>
struct pr{
    f first;
    s second;
    bool is_good;
    pr(f first=0,s second=0):first(first),second(second),is_good(false){}
};

typedef double (*Unary_Function)(double); // unary function
typedef double (*Binary_Function)(double,double); // binary function
typedef double (*Binary_Function_mat)(arma::mat,double); // binary function

template<class RET,class ...Args>
using Mfunction = RET(*)(Args...);

template<class T>
using ConditionFunction = bool(*)(T);
