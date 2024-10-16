//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//[[Rcpp::depends(Rfast)]]

#include "C:\Users\epapadakis\Documents\GitHub\Rfast\inst\include\Rfast\types.hpp"
#include "C:\Users\epapadakis\Documents\GitHub\Rfast\inst\include\Rfast\Set.h"
// using Rfast::Type::Types;
// using Rfast::Type::type;


using namespace std;

#define HASH(key, K)  (3141592653U * (unsigned int)(key) >> (32 - (K)))
#define UTYPEOF(x) ((unsigned)TYPEOF(x))

union uno { double d; unsigned int u[2]; };

template<class T> T* getType(SEXP x){
    if constexpr(is_same<T,SEXP>::value){
        return STRING_PTR(x);
    }else if constexpr(is_same<T,double>::value){
        return REAL(x);
    }else if constexpr(is_same<T,int>::value){
        return INTEGER(x);
    }
    return nullptr;
}

// 
// template <class T, class Func>
// struct hash2
// {
//     int K;
//     Func hash_helper;
//     hash2(int K, Func hash_helper) : K(K), hash_helper(hash_helper) {}
//     std::size_t operator()(const T &key) const
//     {
//         return HASH(hash_helper(key), K);
//     }
// };
// 
// template<class T,class Func> 
// void Table2(SEXP x, const R_xlen_t n, const size_t M, const int K, SEXP &indx, Func hash_helper){
//     size_t ct = 0;
//     hash2<T,Func> hs(K,hash_helper); 
//     std::unordered_set<T,hash2<T,Func>> h(M,hs);
//     T *px = getType<T>(x);
//     for (int i = 0; i < n; ++i)
//     {
//         h.insert(px[i]);
//     }
//     indx = PROTECT(Rf_allocVector(UTYPEOF(x), h.size()));
//     [[maybe_unused]] T *py = nullptr;
//     if constexpr(!is_same<T,SEXP>::value){
//         py = getType<T>(indx);
//     }
//     for (auto v : h)
//     {
//         if constexpr(is_same<T,SEXP>::value){
//             SET_STRING_ELT(indx, ct++, v);
//         }else{
//             py[ct++] = v;
//         }
//     }
//     Rf_copyMostAttrib(x, indx);
//     UNPROTECT(1);
// }
// 
// 
// //[[Rcpp::export]]
// SEXP MyUnique2(SEXP x) {
//     const R_xlen_t n = Rf_xlength(x);
//     const SEXPTYPE tx = UTYPEOF(x);
//     int K;
//     size_t M;
//     SEXP indx = R_NilValue;
//     if (tx == INTSXP || tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
//         if(n >= 1073741824) {
//             stop("Length of 'x' is too large. (Long vector not supported yet)"); // # nocov
//         }
//         const size_t n2 = 2U * (size_t) n;
//         M = 256;
//         K = 8;
//         while (M < n2) {
//             M *= 2;
//             K++;
//         }
//     } else if (tx == LGLSXP) {
//         M = 4;
//         K = 2;
//     } else {
//         stop("Type %s is not supported.", Rf_type2char(tx));
//     }
//     switch(tx){
//     case INTSXP:
//     {
//         Table2<int>(x, n, M, K, indx, [&](int a){ return a; });
//         break;
//     }
//     case REALSXP:
//     {
//         union uno tpv;
//         Table2<double>(x, n, M, K, indx, [&](double a){
//             tpv.d = a;
//             return tpv.u[0] + tpv.u[1];
//         });
//         break;
//     }
//     case STRSXP:
//     {
//         Table2<SEXP>(x, n, M, K, indx, [&](SEXP a){ return (intptr_t)a & 0xffffffff; });
//         break;
//     }
//     }
//     return indx;
// }

template<class T,class Func> 
void group_sum(SEXP x, SEXP ina, const R_xlen_t n, const size_t M, const int K, SEXP &indx, Func hash_helper){
    size_t count = 0, ct = 0, id = 0;
    vector<double> res(n);
    vector<int> pans(n), h(M);
    T *pina = getType<T>(ina);
    double *px = getType<double>(x);
    for (int i = 0; i < n; ++i)
    {
        bool ok = true;
        id = HASH(hash_helper(pina[i]), K);
        while (h[id])
        {
            if (pina[h[id] - 1] == pina[i])
            {
                ok = false;
                break;
            }
            ++id;
            id %= M;
        }
        if(ok){
            ++count;
        }
        if(h[id] == 0){
            h[id] = i + 1;
            pans[i]=i + 1;
        }
        res[h[id]] += px[i];
    }
    vector<int> inds(count);
    for (R_xlen_t i = 0; ct < count; ++i)
    {
        if (pans[i])
        {
            inds[ct++] = pans[i];
        }
    }
    sort(inds.begin(),inds.end(), [&](const int &i,const int &j){
        return pina[i-1] < pina[j-1];
    });
    indx = PROTECT(Rf_allocVector(REALSXP, count));
    double *py = getType<double>(indx);
    for (size_t i = 0; i < count; ++i)
    {
        py[i] = res[inds[i]];
    }
    Rf_copyMostAttrib(x, indx);
    UNPROTECT(1);
}

//[[Rcpp::export]]
SEXP group(SEXP x, SEXP ina) {
    
    // Group<double> g(x,ina,[&](double &a, double &b){ return a+b;})
    
    const R_xlen_t n = Rf_xlength(x);
    const SEXPTYPE tina = UTYPEOF(ina);
    int K;
    size_t M;
    SEXP indx = R_NilValue;
    if (tina == INTSXP || tina == STRSXP || tina == REALSXP || tina == CPLXSXP ) {
        if(n >= 1073741824) {
            stop("Length of 'x' is too large. (Long vector not supported yet)"); // # nocov
        }
        const size_t n2 = 2U * (size_t) n;
        M = 256;
        K = 8;
        while (M < n2) {
            M *= 2;
            K++;
        }
    } else if (tina == LGLSXP) {
        M = 4;
        K = 2;
    } else {
        stop("Type %s is not supported.", Rf_type2char(tina));
    }
    switch(tina){
    case INTSXP:
    {
        group_sum<int>(x, ina, n, M, K, indx, [&](int a){ return a; });
        break;
    }
    case REALSXP:
    {
        union uno tpv;
        group_sum<double>(x, ina, n, M, K, indx, [&](double a){
            tpv.d = a;
            return tpv.u[0] + tpv.u[1];
        });
        break;
    }
    case STRSXP:
    {
        group_sum<SEXP>(x, ina, n, M, K, indx, [&](SEXP a){ return (intptr_t)a & 0xffffffff; });
        break;
    }
    }
    return indx;
}


template<class T, class I> void group_s(SEXP x, SEXP ina, SEXP &indx, const bool sorted){
    auto func = [&](T &a, T &b){ return a+b;};
    // T xx(x);
    Group<T, I, decltype(func)> s(x, ina,func);
    indx = PROTECT(Rf_allocVector(UTYPEOF(x), s.size()));
    s.values(indx, sorted);
    Rf_copyMostAttrib(x, indx);
    UNPROTECT(1);
}

//[[Rcpp::export]]
SEXP group2(SEXP x, SEXP ina, const bool sorted = true) {
    const Types tina = type<SEXP>(ina);
    SEXP indx = R_NilValue;
    switch(tina){
    case Types::INT:
    {
        group_s<int,int>(x, ina, indx, sorted);
        break;
    }
    case Types::REAL:
    {
        group_s<double,double>(x, ina, indx, sorted);
        break;
    }
    // case Types::STRING:
    // {
    //     group_s<NumericVector,SEXP>(x, ina, indx, sorted);
    //     break;
    // }
    default:
        break;
    }
    return indx;
}


template<class T,class Func> 
void Sort_h(SEXP x, const R_xlen_t n, const size_t M, const int K, SEXP &indx, Func hash_helper){
    size_t count = 0, ct = 0, id = 0;
    vector<int> pans(n), h(M), popul(n);
    T *px = getType<T>(x);
    for (int i = 0; i < n; ++i)
    {
        bool ok = true;
        id = HASH(hash_helper(px[i]), K);
        while (h[id])
        {
            if (px[h[id] - 1] == px[i])
            {
                ok = false;
                break;
            }
            ++id;
            id %= M;
        }
        if(ok){
            ++count;
        }
        if(h[id] == 0){
            h[id] = i + 1;
            pans[i]=i+1;
        }
        ++popul[h[id]-1];
    }
    
    //if(count < 0.7 * n) failed
    vector<int> inds(count);
    for (R_xlen_t i = 0; ct < count; ++i)
    {
        if (pans[i])
        {
            inds[ct++] = pans[i]-1;
        }
    }
    sort(inds.begin(),inds.end(), [&](const auto &i,const auto &j){
        return strcmp(CHAR(STRING_ELT(x, i)),CHAR(STRING_ELT(x, j))) < 0; 
    });
    
    indx = PROTECT(Rf_allocVector(UTYPEOF(x), n));
    [[maybe_unused]] T *py = nullptr;
    if constexpr(!is_same<T,SEXP>::value){
        py = getType<T>(indx);
    }
    ct=0;
    for (size_t i = 0; i < count; ++i)
    {
        const auto& size = popul[inds[i]];
        auto& pxx = px[inds[i]];
        // Rcout<<pxx<<" -> "<<v.second<<"\n";
        for(int j = 0; j < size; ++j){
            if constexpr(is_same<T,SEXP>::value){
                SET_STRING_ELT(indx, ct++, pxx);
            }else{
                py[ct++] = pxx;
            }
        }
    }
    Rf_copyMostAttrib(x, indx);
    UNPROTECT(1);
}

//[[Rcpp::export]]
SEXP Sort(SEXP x) {
    const R_xlen_t n = Rf_xlength(x);
    const SEXPTYPE tx = UTYPEOF(x);
    int K;
    size_t M;
    SEXP indx = R_NilValue;
    if (tx == INTSXP || tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
        if(n >= 1073741824) {
            stop("Length of 'x' is too large. (Long vector not supported yet)"); // # nocov
        }
        const size_t n2 = 2U * (size_t) n;
        M = 256;
        K = 8;
        while (M < n2) {
            M *= 2;
            K++;
        }
    } else if (tx == LGLSXP) {
        M = 4;
        K = 2;
    } else {
        stop("Type %s is not supported.", Rf_type2char(tx));
    }
    switch(tx){
    case INTSXP:
    {
        Sort_h<int>(x, n, M, K, indx, [&](int a){ return a; });
        break;
    }
    case REALSXP:
    {
        union uno tpv;
        Sort_h<double>(x, n, M, K, indx, [&](double a){
            tpv.d = a;
            return tpv.u[0] + tpv.u[1];
        });
        break;
    }
    case STRSXP:
    {
        Sort_h<SEXP>(x, n, M, K, indx, [&](SEXP a){ return (intptr_t)a & 0xffffffff; });
        break;
    }
    }
    return indx;
}

template<class T,class Func> 
void freqMax_h(SEXP x, const R_xlen_t n, const size_t M, const int K, SEXP &indx, Func hash_helper){
    size_t count = 0, id = 0;
    vector<int> pans(n), h(M), popul(n);
    T *px = getType<T>(x);
    for (int i = 0; i < n; ++i)
    {
        bool ok = true;
        id = HASH(hash_helper(px[i]), K);
        while (h[id])
        {
            if (px[h[id] - 1] == px[i])
            {
                ok = false;
                break;
            }
            ++id;
            id %= M;
        }
        if(ok){
            ++count;
        }
        if(h[id] == 0){
            h[id] = i + 1;
            pans[i]=i+1;
        }
        ++popul[h[id]-1];
    }
    // auto maxv = max_element(popul.begin(),popul.end());
    // 
    // indx = PROTECT(Rf_allocVector(UTYPEOF(x), 1));
    // if constexpr(is_same<T,SEXP>::value){
    //     SET_STRING_ELT(indx, 0, px[pans[maxv-popul.begin()]-1]);
    // }else{
    //     T *v = getType<T>(indx);
    //     v[0] = px[pans[maxv-popul.begin()]-1];
    // }
    // Rf_copyMostAttrib(x, indx);
    // UNPROTECT(1);
}

//[[Rcpp::export]]
SEXP freqMax(SEXP x) {
    const R_xlen_t n = Rf_xlength(x);
    const SEXPTYPE tx = UTYPEOF(x);
    int K;
    size_t M;
    SEXP indx = R_NilValue;
    if (tx == INTSXP || tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
        if(n >= 1073741824) {
            stop("Length of 'x' is too large. (Long vector not supported yet)"); // # nocov
        }
        const size_t n2 = 2U * (size_t) n;
        M = 256;
        K = 8;
        while (M < n2) {
            M *= 2;
            K++;
        }
    } else if (tx == LGLSXP) {
        M = 4;
        K = 2;
    } else {
        stop("Type %s is not supported.", Rf_type2char(tx));
    }
    switch(tx){
    case INTSXP:
    {
        freqMax_h<int>(x, n, M, K, indx, [&](int a){ return a; });
        break;
    }
    case REALSXP:
    {
        union uno tpv;
        freqMax_h<double>(x, n, M, K, indx, [&](double a){
            tpv.d = a;
            return tpv.u[0] + tpv.u[1];
        });
        break;
    }
    case STRSXP:
    {
        freqMax_h<SEXP>(x, n, M, K, indx, [&](SEXP a){ return (intptr_t)a & 0xffffffff; });
        break;
    }
    }
    return indx;
}


template<class T,class Func> 
void Table_h(SEXP x, const R_xlen_t n, const size_t M, const int K, SEXP &indx, Func hash_helper){
    size_t count = 0, ct = 0, id = 0;
    vector<int> pans(n), h(M), popul(n);
    T *px = getType<T>(x);
    for (R_xlen_t i = 0; i < n; ++i)
    {
        bool ok = true;
        id = HASH(hash_helper(px[i]), K);
        while (h[id])
        {
            if (px[h[id] - 1] == px[i])
            {
                ok = false;
                break;
            }
            ++id;
            id %= M;
        }
        if(ok){
            ++count;
        }
        if(h[id] == 0){
            h[id] = i + 1;
            pans[i]=i+1;
        }
        ++popul[h[id]-1];
    }
    
    indx = PROTECT(Rf_allocVector(INTSXP, count));
    int *py = getType<int>(indx);
    for (R_xlen_t i = 0; ct < count; ++i)
    {
        if(pans[i]){
            py[ct++] = popul[pans[i]-1];
        }
    }
    Rf_copyMostAttrib(x, indx);
    UNPROTECT(1);
}

//[[Rcpp::export]]
SEXP Tablert(SEXP x) {
    const R_xlen_t n = Rf_xlength(x);
    const SEXPTYPE tx = UTYPEOF(x);
    int K;
    size_t M;
    SEXP indx = R_NilValue;
    if (tx == INTSXP || tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
        if(n >= 1073741824) {
            stop("Length of 'x' is too large. (Long vector not supported yet)"); // # nocov
        }
        const size_t n2 = 2U * (size_t) n;
        M = 256;
        K = 8;
        while (M < n2) {
            M *= 2;
            K++;
        }
    } else if (tx == LGLSXP) {
        M = 4;
        K = 2;
    } else {
        stop("Type %s is not supported.", Rf_type2char(tx));
    }
    switch(tx){
    case INTSXP:
    {
        Table_h<int>(x, n, M, K, indx, [&](int a){ return a; });
        break;
    }
    case REALSXP:
    {
        union uno tpv;
        Table_h<double>(x, n, M, K, indx, [&](double a){
            tpv.d = a;
            return tpv.u[0] + tpv.u[1];
        });
        break;
    }
    case STRSXP:
    {
        Table_h<SEXP>(x, n, M, K, indx, [&](SEXP a){ return (intptr_t)a & 0xffffffff; });
        break;
    }
    }
    return indx;
}

template<class T>
void parallelQuickSort(T& arr, int low, int high, int maxDepth) {
    if (low < high) {
        if (maxDepth <= 0) {
            std::sort(arr.begin() + low, arr.begin() + high);
        } else {
            int pivotIndex = low + (high - low) / 2;
            int pivotValue = arr[pivotIndex];
            std::swap(arr[pivotIndex], arr[high - 1]);
            int storeIndex = low;
            for (int i = low; i < high - 1; i++) {
                if (arr[i] < pivotValue) {
                    std::swap(arr[storeIndex], arr[i]);
                    storeIndex++;
                }
            }
            std::swap(arr[storeIndex], arr[high - 1]);
#pragma omp task shared(arr) if (high - low > 1000)
            parallelQuickSort(arr, low, storeIndex - 1, maxDepth - 1);
#pragma omp task shared(arr) if (high - low > 1000)
            parallelQuickSort(arr, storeIndex + 1, high, maxDepth - 1);
#pragma omp taskwait
        }
    }
}

template<class T>
void parallelIntrosort(T& arr, int maxDepth = 1) {
#pragma omp parallel
#pragma omp single
    parallelQuickSort(arr, 0, arr.size(), maxDepth);
}

//[[Rcpp::export]]
SEXP ssort(SEXP x, int maxDepth = 0){
    SEXP xx = Rf_duplicate(x);
    const SEXPTYPE tx = UTYPEOF(x);
    const int n = Rf_length(x);
    switch(tx){
        case INTSXP:
        {
            arma::icolvec X(INTEGER(xx),n,false);
            parallelIntrosort(X,maxDepth * std::log2(n));
            break;
        }
        case REALSXP:
        {
            arma::colvec X(REAL(xx),n,false);
            parallelIntrosort(X,maxDepth * std::log2(n));
            break;
        }
        case STRSXP:
        {
            break;
        }
    }
    return xx;
}
