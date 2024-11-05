#pragma once

#ifdef WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#endif
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rversion.h>
#if !defined(R_VERSION) || R_VERSION < R_Version(3, 5, 0)
#define USE_RINTERNALS
#define DATAPTR_RO(x) ((const void *)DATAPTR(x))
#endif
#include <Rinternals.h>

#include "types.hpp"
#include <vector>

using namespace Rfast::Type;
using std::is_same;
using std::vector;

#define N_ISNAN(x, y) (!ISNAN(x) && !ISNAN(y))
#define B_IsNA(x, y) (R_IsNA(x) && R_IsNA(y))
#define B_IsNaN(x, y) (R_IsNaN(x) && R_IsNaN(y))
#define B_ISNAN(x, y) (ISNAN(x) && ISNAN(y))
#define C_IsNA(x) (R_IsNA(x.r) || R_IsNA(x.i))
#define C_IsNaN(x) (R_IsNaN(x.r) || R_IsNaN(x.i))
#define C_ISNAN(x, y) (B_ISNAN(x, y) || (N_ISNAN(x, y) && x == y))

template <class T>
class Hash
{
protected:
    size_t M = 256, un_len = 0;
    int K = 8;
    unsigned int HASH(unsigned int key)
    {
        return (3141592653U * (key) >> (32 - (K)));
    }

    void initMemory(SEXP x, unsigned int n)
    {
        const Types tx = type<T>(x);
        if (tx == Types::INT || tx == Types::STRING || tx == Types::REAL || tx == Types::COMPLEX)
        { //
            if (n >= 1073741824)
            {
                stop("Length of 'x' is too large. (Long vector not supported yet)"); // # nocov
            }
            const size_t n2 = 2U * (size_t)n;
            while (M < n2)
            {
                M <<= 1; // M *=2
                K++;
            }
        }
        else if (tx == Types::LOGICAL)
        {
            M = 4;
            K = 2;
        }
    }
};

template <class T>
struct HashBase
{
    inline T operator()(T &a) { return a; }
};

template <>
struct HashBase<SEXP>
{
    inline intptr_t operator()(SEXP &a) { return (intptr_t)a & 0xffffffff; }
};

template <>
struct HashBase<double>
{
    union uno
    {
        double d;
        unsigned int u[2];
    } tpv;
    inline unsigned int operator()(double &a)
    {
        tpv.d = R_IsNA(a) ? Rfast::NA<double>::value() : (R_IsNaN(a) ? R_NaN : a);
        return tpv.u[0] + tpv.u[1];
    }
};

template <>
struct HashBase<Rcomplex>
{
    union uno
    {
        double d;
        unsigned int u[2];
    } tpv;
    Rcomplex tmp;
    inline unsigned int operator()(Rcomplex &a)
    {
        tmp.r = (a.r == 0.0) ? 0.0 : a.r;
        tmp.i = (a.i == 0.0) ? 0.0 : a.i;
        if (C_IsNA(tmp))
        {
            tmp.r = tmp.i = Rfast::NA<double>::value();
        }
        else if (C_IsNaN(tmp))
        {
            tmp.r = tmp.i = R_NaN;
        }
        tpv.d = tmp.r;
        unsigned int u = tpv.u[0] ^ tpv.u[1];
        tpv.d = tmp.i;
        u ^= tpv.u[0] ^ tpv.u[1];
        return u;
    }

    friend bool operator==(Rcomplex &x, Rcomplex &y)
    {
        return (N_ISNAN(x.r, x.i) && N_ISNAN(y.r, y.i)) ? (x.r == y.r && x.i == y.i) : (C_IsNA(x) ? C_IsNA(y) : (C_IsNA(y) ? 0 : (C_ISNAN(x.r, y.r) && C_ISNAN(x.i, y.i))));
    }
};

template <class T, class Hash_Helper = HashBase<T>>
class Set : protected Hash<T>
{
    T *data;
    size_t n, un_len = 0;
    Hash_Helper hash_helper;
    vector<int> pans, h;

    void insertAll(const bool fromLast = false)
    {
        if (fromLast)
        {
            for (size_t i = n; i > 0; --i)
            {
                insert(data[i - 1], i - 1);
            }
        }
        else
        {
            for (size_t i = 0; i < n; ++i)
            {
                insert(data[i], i);
            }
        }
    }

public:
    template <class S>
    Set(S x, const bool fromLast = false, Hash_Helper hash_helper = HashBase<T>()) : hash_helper(hash_helper)
    {
        if constexpr (is_same<S, SEXP>::value)
        {
            this->data = Rfast::asPtr<T>(x);
            this->n = Rf_length(x);
        }
        else
        {
            this->data = x.data();
            this->n = x.size();
        }
        this->initMemory(x, n);
        pans = vector<int>(n);
        h = vector<int>(this->M);
        insertAll(fromLast);
    }

    inline void insert(T &d, int i = 0)
    {
        bool ok = true;
        size_t id = this->HASH(hash_helper(d));
        while (h[id])
        {
            if (data[h[id] - 1] == d)
            {
                ok = false;
                break;
            }
            ++id;
            id %= this->M;
        }
        if (ok)
        {
            h[id] = i + 1;
            ++pans[i];
            ++un_len;
        }
    }

    size_t size() const
    {
        return un_len;
    }

    inline void values(SEXP &indx)
    {
        [[maybe_unused]] T *py = nullptr;
        if constexpr (!is_same<T, SEXP>::value)
        {
            py = Rfast::asPtr<T>(indx);
        }
        size_t ct = 0;
        for (size_t i = 0; ct < un_len; ++i)
        {
            if (pans[i])
            {
                if constexpr (is_same<T, SEXP>::value)
                {
                    SET_STRING_ELT(indx, ct++, data[i]);
                }
                else
                {
                    py[ct++] = data[i];
                }
            }
        }
    }
};

template <class T, class I, class Function, class Hash_Helper = HashBase<I>>
class Group : protected Hash<I>
{
    T *data;
    I *ina;
    size_t n, un_len = 0;
    Function func;
    Hash_Helper hash_helper;
    vector<int> pans, h;
    vector<T> res;
    using result_type = T; // typename std::remove_reference<typename T::value_type>::type;

public:
    Group(SEXP x, SEXP ina, Function func, T init_val = 0, Hash_Helper hash_helper = Hash_Helper()) : func(func), hash_helper(hash_helper)
    {
        this->data = Rfast::asPtr<T>(x);
        this->ina = Rfast::asPtr<I>(ina);
        this->n = Rf_length(x);
        this->initMemory(ina, n);
        pans = vector<int>(n);
        h = vector<int>(this->M);
        res = vector<T>(n, init_val);
        for (size_t i = 0; i < n; ++i)
        {
            insert(this->ina[i], i);
        }
    }

    inline void insert(I &d, int i = 0)
    {
        bool ok = true;
        size_t id = this->HASH(hash_helper(d));
        while (h[id])
        {
            if (ina[h[id] - 1] == d)
            {
                ok = false;
                break;
            }
            ++id;
            id %= this->M;
        }
        if (ok)
        {
            ++un_len;
        }
        if (h[id] == 0)
        {
            h[id] = i + 1;
            pans[i] = i + 1;
        }
        res[h[id]] = func(res[h[id]], data[i]);
    }

    size_t size() const
    {
        return un_len;
    }

    inline void values(SEXP &indx, const bool sorted = true)
    {
        size_t ct = 0;
        vector<int> inds(un_len);
        for (size_t i = 0; ct < un_len; ++i)
        {
            if (pans[i])
            {
                inds[ct++] = pans[i];
            }
        }
        if (sorted)
        {
            sort(inds.begin(), inds.end(), [&](const int &i, const int &j)
                 { return ina[i - 1] < ina[j - 1]; });
        }
        result_type *py = Rfast::asPtr<result_type>(indx);
        for (size_t i = 0; i < un_len; ++i)
        {
            py[i] = res[inds[i]];
        }
    }
};

    
template <class T, class I, class Hash_Helper = HashBase<I>>
class GroupBucket : protected Hash<I>
{
public:
    using Bucket = vector<T>;

private:
    T *data;
    I *ina;
    size_t n, un_len = 0;
    Hash_Helper hash_helper;
    vector<int> pans, h;
    vector<Bucket> res;
    using result_type = T; // typename std::remove_reference<typename T::value_type>::type;

public:

    GroupBucket(SEXP x, SEXP ina, T init_val = 0, Hash_Helper hash_helper = Hash_Helper()) : hash_helper(hash_helper)
    {
        this->data = Rfast::asPtr<T>(x);
        this->ina = Rfast::asPtr<I>(ina);
        this->n = Rf_length(x);
        this->initMemory(ina, n);
        pans = vector<int>(n);
        h = vector<int>(this->M);
        res = vector<Bucket>(n, Bucket(0,init_val));
        for (size_t i = 0; i < n; ++i)
        {
            insert(this->ina[i], i);
        }
    }

    inline void insert(I &d, int i = 0)
    {
        bool ok = true;
        size_t id = this->HASH(hash_helper(d));
        while (h[id])
        {
            if (ina[h[id] - 1] == d)
            {
                ok = false;
                break;
            }
            ++id;
            id %= this->M;
        }
        if (ok)
        {
            ++un_len;
        }
        if (h[id] == 0)
        {
            h[id] = i + 1;
            pans[i] = i + 1;
        }
        res[h[id]].emplace_back(data[i]);
    }

    size_t size() const
    {
        return un_len;
    }

    template<class Function>
    inline void values(SEXP &indx, const bool sorted, Function func)
    {
        size_t ct = 0;
        vector<int> inds(un_len);
        for (size_t i = 0; ct < un_len; ++i)
        {
            if (pans[i])
            {
                inds[ct++] = pans[i];
            }
        }
        if (sorted)
        {
            sort(inds.begin(), inds.end(), [&](const int &i, const int &j)
                 { return ina[i - 1] < ina[j - 1]; });
        }
        result_type *py = Rfast::asPtr<result_type>(indx);
        for (size_t i = 0; i < un_len; ++i)
        {
            py[i] = func(res[inds[i]]);
        }
    }
};
