// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <R.h>
#include <Rinternals.h>


//#include <Rinlinedfuns.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

#ifndef TEMPLATES_H
#define TEMPLATES_H

//[[Rcpp::plugins(cpp11)]]

struct NA {
  template<class T> static T value(){
    return is_same<T,int>::value ? NA_INTEGER : NA_REAL;
  }
};

template<typename f,typename s>
struct pr{
  f first;
  s second;
  bool is_good;
  pr(f first=0,s second=0):first(first),second(second),is_good(false){}
};

typedef double (*Unary_Function)(double); // unary function
typedef double (*Binary_Function)(double,double); // binary function
typedef double (*Binary_Function_mat)(mat,double); // binary function

template<class RET,class ...Args>
using Mfunction = RET(*)(Args...);

template<class T>
using ConditionFunction = bool(*)(T);

/*
 * F: unary function
 * T1,T2: arguent class
 *
 * fill memory that starts from "startf" applying function F from "start" to "end"
*/
template<Unary_Function F,typename T1,typename T2>
void fill_with(T1 start,T1 end,T2 startf){
  for(;start!=end;++start,++startf)
    *startf=F(*start);
}


/*
 * F: unary function
 * T1,T2: arguent class
 *
 * fill memory that starts from "startf" applying function F to "x"
*/
template<Unary_Function F,typename T1,typename T2>
void fill_with(T1 x,T2 startf){
  for(typename T1::iterator start=x.begin();start!=x.end();++start)
    *startf=F(*start);
}

/*
 * F: unary function
 * T1,T2: arguent class
 *
 * fill memory that starts from "startf" applying function F to "x"
*/
template<Unary_Function F,ConditionFunction<double> COND,typename T1,typename T2>
void fill_with_condition(T1 start,T1 end,T2 startf){
  for(;start!=end;++start,++startf)
    if(COND(*start))
      *startf=F(*start);
}


/*
 * T: argument and return class
 *
 * applying function F to "x" and return "x".
*/
template<Unary_Function F,typename T>
T foreach(T x){
  for(typename T::iterator start=x.begin();start!=x.end();++start)
      *start=F(*start);
    return x;
}


/*
 * T1,T2: argument class
 *
 * count the frequency of "value" in "x"
*/
template<typename T1,typename T2 >
int count_value_helper(T1 x,T2 value){
  int s=0;
  for(typename T1::iterator start=x.begin();start!=x.end();++start)
    if(*start==value)
      ++s;
  return s;
}


/*
 * T: argument class
 *
 * count the frequency of "NAs" in "x"
*/
template<typename T>
int count_NAs(T x){
  int s=0;
  for(typename T::iterator start=x.begin();start!=x.end();++start)
    if(R_IsNA(*start))
      ++s;
  return s;
}


/*
 * T: argument class
*/
template<typename T>
double med_helper(typename T::iterator start,typename T::iterator last){
  double F;
  const int sz=last-start,middle=sz/2-1;
  if(sz%2==0){
    std::nth_element(start,start+middle,last);
    F=(start[middle]+*(std::min_element(start+middle+1,last)))/2.0;
  }else{
    std::nth_element(start,start+middle+1,last);
    F=start[middle+1];
  }
  return F;
}


/*
 * Ret: return class
 * T: argument class
*/
template<typename Ret,typename T>
Ret Order(T x,const bool stable,const bool descend,const int init_v){
  Ret ind(x.size());
  iota(ind.begin(),ind.end(),init_v);
  if(descend){
    auto descend_func = [&](int i,int j){return x[i-init_v]>x[j-init_v];};
    stable ? stable_sort(ind.begin(),ind.end(),descend_func) : sort(ind.begin(),ind.end(),descend_func);
  }else{
    auto func = [&](int i,int j){return x[i-init_v]<x[j-init_v];};
    stable ? stable_sort(ind.begin(),ind.end(),func) : sort(ind.begin(),ind.end(),func);
  }
  return ind;
}

/*
 * T: argument class (not classes)
*/
template<typename T>
void min_max(T *start,T *end,T &mn, T &mx){
  T xxx;
  mn=mx=*start;
  start++;
  for(;start!=end;++start){
    xxx=*start;
    if(xxx>mx)
      mx=xxx;
    else if(xxx<mn)
      mn=xxx;
  }
}

/*
 * T: argument class (not classes)
*/
template<typename T>
void maximum(T *start,T *end,T &mx){
  double xxx;
  mx=*start;
  start++;
  for(;start!=end;++start){
    xxx=*start;
    if(xxx>mx)
      mx=xxx;
  }
}


/*
 * T: argument class (not classes)
*/
template<typename T>
void minimum(T *start,T *end,T &mn){
  double xxx;
  mn=*start;
  start++;
  for(;start!=end;++start){
    xxx=*start;
    if(xxx<mn)
      mn=xxx;
  }
}


/*
 * Ret: return class
 * T: argument class
*/
template<typename Ret,typename T>
Ret Order_rank(T& x,const bool descend,const bool stable,const int n,const int k){
  Ret ind(x.size()-k);
  iota(ind.begin(),ind.end(),0);
  if(descend){
    auto descend_func = [&](int i,int j){return x[i]>x[j];};
    stable ? stable_sort(ind.begin(),ind.end()-n,descend_func) : sort(ind.begin(),ind.end()-n,descend_func);
  }else{
    auto func = [&](int i,int j){return x[i]<x[j];};
    stable ? stable_sort(ind.begin(),ind.end()-n,func) : sort(ind.begin(),ind.end()-n,func);
  }
  return ind;
}

/*
 * F: unary function
 * T1,T2: arguent class
 *
 * fill memory that starts from "startf" applying function F from "start" to "end"
*/
template<ConditionFunction<double> COND,class T>
double sum_with_condition(T x){
   double a=0.0;
  for(typename T::iterator start=x.begin();start!=x.end();++start)
    if(COND(*start))
      a+=*start;
  return a;
}

/*
 * F: a unary function 
 * T: argument class
 *
 * applying function F to "x" and sum the elements of "x"
*/
template<Unary_Function F,typename T>
double sum_with(T x){
  double a=0;
  for(typename T::iterator start=x.begin();start!=x.end();++start)
    a+=F(*start);
  return a;
}

/*
 * T: argument class
 *
 * applying function F to "x" and sum the elements of "x"
*/
template<typename T,Mfunction<T,T> F>
double sum_with(T *start,T *end){
  double a=0;
  for(;start!=end;++start)
    a+=F(*start);
  return a;
}

/*
 * F: a binary function 
 * T: argument class
 *
 * applying function F to "x" and sum the elements of "x"
*/
template<Binary_Function F,typename T>
double sum_with(T x,const double p){
  double a=0;
  for(typename T::iterator start=x.begin();start!=x.end();++start)
    a+=F(*start,p);
  return a;
}

/*
 * F: a binary function 
 * T: argument class
 *
 * applying function F to "x" and sum the elements of "x"
*/
template<Binary_Function F,class T1,class T2>
double sum_with(T1 x,T2& y){
  double a=0;
  typename T1::iterator startx=x.begin();
  typename T2::iterator starty=y.begin();
  for(;startx!=x.end();++startx,++starty)
    a+=F(*startx,*starty);
  return a;
}

/*
 * F1: a binary function
 * F2: a binary function
 * T1: argument class
 * T2: argument class
 *
 * applying function F1 to "x","y" and applying F2 to the result
*/
template<class T1,class T2>
double Apply(T1 x,T2& y,Binary_Function F1,Binary_Function F2){
  double a=0;
  typename T1::iterator startx=x.begin();
  typename T2::iterator starty=y.begin();
  for(;startx!=x.end();++startx,++starty)
    a=F2(a,F1(*startx,*starty));
  return a;
}

/*
 * F1: a binary function
 * F2: a binary function
 * T1: argument class
 * T2: argument class
 *
 * applying function F1 to "x","y" and applying F2 to the result
*/
template<class T1,class T2,Binary_Function F1,Binary_Function F2>
double Apply(T1 x,T2& y){
  double a=0;
  typename T1::iterator startx=x.begin();
  typename T2::iterator starty=y.begin();
  for(;startx!=x.end();++startx,++starty)
    a=F2(a,F1(*startx,*starty));
  return a;
}

/*
 * F1: a binary function
 * F2: a binary function
 * T1: argument class
 * T2: argument class
 *
 * applying function F1 to "x","y" and applying F2 to the result
*/
template<class T1,class T2,Binary_Function F1,Binary_Function F2>
double Apply(T1 *x,T1 *endx,T2 *y){
  double a=0;
  for(;x!=endx;++x,++y)
    a=F2(a,F1(*x,*y));
  return a;
}

/*
 * F1: a binary function
 * F2: a binary function
 * T1: argument class
 * T2: argument class
 *
 * applying function F1 to "x","y" and applying F2 to the result
*/
template<class T1,class T2,Binary_Function F1,Unary_Function F2,Binary_Function F3>
double Apply(T1 x,T2& y){
  double a=0;
  typename T1::iterator startx=x.begin();
  typename T2::iterator starty=y.begin();
  for(;startx!=x.end();++startx,++starty)
    a=F3(a,F2(F1(*startx,*starty)));
  return a;
}

/*
 * F1: a binary function
 * F2: a binary function
 * T1: argument class
 * T2: argument class
 *
 * applying function F1 to "x","y" and applying F2 to the result
*/
template<class T1,class T2>
double Apply(T1 x,T2& y,Binary_Function F1,Unary_Function F2,Binary_Function F3){
  double a=0;
  typename T1::iterator startx=x.begin();
  typename T2::iterator starty=y.begin();
  for(;startx!=x.end();++startx,++starty)
    a=F3(a,F2(F1(*startx,*starty)));
  return a;
}

/*
 * F1: a binary function
 * F2: a binary function
 * T1: argument class
 * T2: argument class
 *
 * applying function F1 to "x","y" and applying F2 to the result
*/
template<class T,Unary_Function F1,Binary_Function F2>
double Apply(T x){
  double a=0;
  typename T::iterator start=x.begin();
  for(;start!=x.end();++start)
    a=F2(a,F1(*start));
  return a;
}

/*
 * T: return and argument type (Rcpp and simple types. Not armadillo)
*/
template<typename T>
T square2(T x){
  return x*x;
}

/*
 * Ret: return class
 * T: argument class
*/
template<typename Ret,typename T>
Ret Tabulate(T x,int &nroww){
  Ret f(nroww);
  std::fill(f.begin(),f.end(),0);
  typename Ret::iterator F=f.begin();
  typename T::iterator a=x.begin();
  for(;a!=x.end();++a){
      F[*a-1]++;
  }
  return f;
}

/*
 * Ret: return class
 * T: argument class
*/
template<typename Ret>
Ret Tabulate(int* start,int* end,int& nroww){
  Ret f(nroww);
  std::fill(f.begin(),f.end(),0);
  typename Ret::iterator F=f.begin();
  for(;start!=end;++start){
    F[*start-1]++;
  }
  return f;
}

/*
 * Ret: return class
 * T1: argument class
 * T2: argument class
*/
/*
 * Ret: return class
 * T1: argument class
 * T2: argument class
*/
template<typename Ret,typename T1,typename T2>
Ret group_sum_helper(T1 x,T2 key,int *minn=nullptr,int *maxx=nullptr){
  int mn,mx;
  const bool is_mn_null=(minn==nullptr),is_mx_null=(maxx==nullptr);
  if(is_mx_null && is_mn_null){
    min_max<int>(key.begin(),key.end(),mn,mx);
  }else if(is_mx_null){
    mn=*minn;
    maximum<int>(key.begin(),key.end(),mx);
  }else if(is_mn_null){
    mx=*maxx;
    minimum<int>(key.begin(),key.end(),mn);
  }else{
    mn=*minn;
    mx=*maxx;
  }
  typename T2::iterator kk=key.begin();
  vector<double> f(mx-mn+1);
  vector<bool> is_good(mx-mn+1);
  typename T1::iterator xx=x.begin(),rr;
  vector<double>::iterator ff=f.begin();
  vector<bool>::iterator ok;
  int index;
  for(;xx!=x.end();++xx,++kk){
    index = *kk-mn;
    is_good[index]=true;
    f[index] += *xx;
  }
  int number_of_values=0;
  for(auto v : is_good){
    if(v)
      ++number_of_values;
  }
  Ret res(number_of_values);
  for(rr=res.begin(),ok=is_good.begin(),ff=f.begin();ff!=f.end();++ff){
    if(*ok++)
      *rr++=*ff;
  }
  return res;
}


/* Distixos gia tora to T prepei na einai panta mat 
 * dioti mono ta class tou arma exoun n_rows,n_cols
 * Ret: return class
 * T: argument class
*/
template<typename Ret,typename T,typename Vec>
Ret cross_x_y(T x,T y){
  const int ncl=x.n_cols,nrw=x.n_rows,p=y.n_cols;
  Ret f(ncl,p);
  Vec yv(nrw);
  int i,j;
  for(i=0;i<p;++i){
    yv=y.col(i);
    for(j=0;j<ncl;++j){
      f(j,i)=sum(x.col(j)%yv);
    }
  }
  return f;
}

/*
 * Ret: return class
 * T: argument class
*/
template<typename Ret,typename T>
Ret cross_x(T x){
  const int ncl=x.n_cols;
  Ret f(ncl,ncl);
  double a;
  int i,j;
  for(i=0;i<ncl;++i){
    for(j=i;j<ncl;++j){
      a=sum(x.col(j)%x.col(i));
      f(i,j)=a;
      f(j,i)=a;
    }
  }

    return f;
}

/*
 * kanonika to T prepei na einai class tou Rcpp allios h sort_unique tha vgalei error
 * Ret: return class
 * T: argument class
*/
template<typename Ret,typename T>
Ret design_matrix_helper(T x) {
  int i=0;
  const int n=x.size();
  T tmp=sort_unique(x);
  typename T::iterator xx=x.begin(),leksi_bg,leksi_en;
  Ret Final(n,tmp.size());
  std::fill(Final.begin(),Final.end(),0);
  for(leksi_bg=tmp.begin(),leksi_en=tmp.end(),i=0;xx!=x.end();++xx,++i)
    Final(i,lower_bound(leksi_bg,leksi_en,*xx)-leksi_bg)=1;
  return Final;
}

template<typename T>
inline T madd(T x,T y){
  return x+y;
}


template<typename T>
inline T mless(T x,T y){
  return x<y;
}

template<typename T>
inline T mless_eq(T x,T y){
  return x<=y;
}

template<typename T>
inline T mgreater_eq(T x,T y){
  return x>=y;
}

template<typename T>
inline T mgreater(T x,T y){
  return x>y;
}


template<typename T>
inline T mdiv(T x,T y){
  return x/y;
}


template<typename T>
inline T mdiff(T x,T y){
  return x-y;
}


template<typename T>
inline T mmult(T x,T y){
  return x*y;
}

template<typename T>
inline T mequal(T x,T y){
  return x==y;
}

template<typename T>
inline T mmax(T x,T y){
  return std::max(x,y);
}

template<typename T>
inline T mmin(T x,T y){
  return std::min(x,y);
}

template<typename T>
inline T mand(T x,T y){
  return x&&y;
}

template<typename T>
inline T mor(T x,T y){
  return x||y;
}

template<typename RET,typename T1,typename T2>
inline RET madd(T1 x,T2 y){
  return x+y;
}


template<typename RET,typename T1,typename T2>
inline RET mless(T1 x,T2 y){
  return x<y;
}

template<typename RET,typename T1,typename T2>
inline RET mless_eq(T1 x,T2 y){
  return x<=y;
}

template<typename RET,typename T1,typename T2>
inline RET mgreater_eq(T1 x,T2 y){
  return x>=y;
}

template<typename RET,typename T1,typename T2>
inline RET mgreater(T1 x,T2 y){
  return x>y;
}


template<typename RET,typename T1,typename T2>
inline RET mdiv(T1 x,T2 y){
  return x/y;
}


template<typename RET,typename T1,typename T2>
inline RET mdiff(T1 x,T2 y){
  return x-y;
}


template<typename RET,typename T1,typename T2>
inline RET mmult(T1 x,T2 y){
  return x*y;
}

template<typename RET,typename T1,typename T2>
inline RET mequal(T1 x,T2 y){
  return x==y;
}

template<typename RET,typename T1,typename T2>
inline RET mmax(T1 x,T2 y){
  return std::max(x,y);
}

template<typename RET,typename T1,typename T2>
inline RET mmin(T1 x,T2 y){
  return std::min(x,y);
}

template<typename RET,typename T1,typename T2>
inline RET mand(T1 x,T2 y){
  return x&&y;
}

template<typename RET,typename T1,typename T2>
inline RET mor(T1 x,T2 y){
  return x||y;
}


template<typename T, Binary_Function F>
void myoperator(T f[],T &x,T *y,int &len){
  T *ff=f;
  for(int i=0;i<len;++i,++ff,++y)
    *ff=F(x,*y);
}

template<typename T>
void as_integer_h_sorted(vector<T> x,IntegerVector &f,const int init,const T val){
  const int n=x.size();
  int i,j=0,c=init;
  sort(x.begin(),x.end());
  auto v=x[j];
  f[0]=init;
  for(i=1;i<n;++i){
    if(v!=x[i]){
      j=i;
      v=x[j];
      ++c;
    }
    f[i]=c;
  }
}


template<typename T>
void as_integer_h(vector<T> x,IntegerVector &f,const int init,const T val){
  const int n=x.size();
  int i,j=0,c=init;
  vector<int> ind=Order< vector<int>,vector<T> >(x,false,false,0); // diorthoseiii
  x.push_back(val);
  T v=x[ind[j]];
  f[ind[0]]=init;
  for(i=1;i<n;++i){
    if(v!=x[ind[i]]){
      j=i;
      v=x[ind[j]];
      ++c;
    }
    f[ind[i]]=c;
  }
}

template<typename T>
void as_integer_h_with_names(vector<T> x,List &L,const int init,const T val){
  const int n=x.size()+1;
  int i,j=0,c=init;
  vector<int> ind=Order< vector<int>,vector<T> >(x,false,false,0); // diorthoseiii
  x.push_back(val);
  ind.push_back(0);
  vector<T> w;
  T v=x[ind[j]];
  IntegerVector f(n-1);
  f[ind[0]]=init;
  for(i=1;i<n;++i){
    if(v!=x[ind[i]]){
      j=i;
      w.push_back(v);
      v=x[ind[j]];
      ++c;
    }
    f[ind[i]]=c;
  }
  L["w"]=w;
  L["f"]=f;
}

template<typename T>
vector<int> table_use_na(vector<T> x,const int use_na){
  typename vector<T>::iterator new_end = std::remove_if(x.begin(),x.end(),R_IsNA);
  sort(x.begin(),new_end);
  vector<int> n;
  typename vector<T>::iterator a=x.begin(),b=a+1;
  List l;
  int s=1;
  for(;b!=new_end;++b){
  if(*a!=*b){
    n.push_back(s);
    a=b;
    s=1;
  }else
    ++s;
  }
  if(use_na==1){
    n.push_back(new_end-x.begin());
  }
  return n;
}
 
template<typename T>
vector<int> table_simple(vector<T> x){
  sort(x.begin(),x.end());
  x.push_back(T(0));
  vector<int> n;
  typename vector<T>::iterator a=x.begin(),b=a+1;
  int s=1;
  for(;b!=x.end();++b){
    if(*a!=*b){
      n.push_back(s);
      a=b;
      s=1;
    }else
      ++s;
  }
  return n;
}

template<typename T>
void table2_like_r(vector<T> x,vector<T> y,IntegerMatrix &f,const T val){
  //cout<<"table2_like_r\n";
  const int n=x.size();
  int mx_x,mx_y;
  IntegerVector ix(n),iy(n);
  as_integer_h<T>(x,ix,0,val);
  as_integer_h<T>(y,iy,0,val);
  maximum<int>(ix.begin(),ix.end(),mx_x);
  maximum<int>(iy.begin(),iy.end(),mx_y);
  f=IntegerMatrix(mx_x+1,mx_y+1);
  for(int i=0;i<n;++i){
    f(ix[i],iy[i])++; // den xreiazetai -1 efoson kaleitai h as_integer me 0
  }
}


template<typename T>
void table2_like_r_with_names(vector<T> x,vector<T> y,List &L,const T val){
  //cout<<"table2_like_r_with_names"<<endl;
  const int n=x.size();
  int mx_x,mx_y;
  List lx,ly;
  as_integer_h_with_names<T>(x,lx,0,val);
  as_integer_h_with_names<T>(y,ly,0,val);
  IntegerVector ix=lx["f"],iy=ly["f"];
  maximum<int>(ix.begin(),ix.end(),mx_x);
  maximum<int>(iy.begin(),iy.end(),mx_y);
  IntegerMatrix f(mx_x+1,mx_y+1);
  for(int i=0;i<n;++i){
    f(ix[i],iy[i])++; // den xreiazetai -1 efoson kaleitai h as_integer me 0
  }
  L["x"]=lx["w"];
  L["y"]=ly["w"];
  L["f"]=f;
}

template<typename Ret,typename T,typename I>
Ret rank_mean(T x,const bool descend){
  const int n=x.size(),n_1=n+1;
  int i,j=0;
  x.resize(n_1);
  x[n]=0;
  I ind=Order_rank<I,T>(x,descend,false,1,0);
  int k=0,m,times=0;
  Ret f(n);
  double mn=0.0,v=x[ind[j]];
  for(i=1;i<n_1;++i){
    if(v!=x[ind[i]]){
      times=i-j;
      mn=(j+1+i)*0.5; //mn=mean(seq(j+1,i));
      for(k=j,m=0;m<times;++m)
        f[ind[k++]]=mn;
      j=i;
      v=x[ind[j]];
    }
  }
  return f;
}

template<typename Ret,typename T,typename I>
Ret rank_max(T x,const bool descend){
  const int n=x.size(),n_1=n+1;
  int i,j=0;
  x.resize(n_1);
  x[n]=0;
  I ind=Order_rank<I,T>(x,descend,false,1,0);
  int k=0,m,times=0;
  Ret f(n);
  double v=x[ind[j]];
  for(i=1;i<n_1;++i){
    if(v!=x[ind[i]]){
      times=i-j;
      for(k=j,m=0;m<times;++m)
        f[ind[k++]]=i;
      j=i;
      v=x[ind[j]];
    }
  }
  return f;
}

template<typename Ret,typename T,typename I>
Ret rank_min(T x,const bool descend){
  const int n=x.size();
  int i,j=0;
  I ind=Order_rank<I,T>(x,descend,false,0,1);
  Ret f(n);
  double v=x[ind[j]];
  f[ind[0]]=1;
  for(i=1;i<n;++i){
    if(v!=x[ind[i]]){
      j=i;
      v=x[ind[j]];
    }
    f[ind[i]]=j+1;
  }
  return f;
}

template<typename Ret,typename T,typename I>
Ret rank_first(T x,const bool descend,const bool stable){
  const int n=x.n_elem;
  I ind=Order_rank<I,T>(x,descend,stable,0,1);
  Ret f(n);
  for(int i=0;i<n;++i){
    f[ind[i]]=i+1;
  }
  return f;
}

template<Binary_Function oper,Binary_Function func>
NumericVector eachcol_apply_helper(NumericMatrix& x,NumericVector& y,SEXP ind){
  const bool is_ind_null = Rf_isNull(ind);
  const int n = is_ind_null ? x.ncol() : LENGTH(ind);
  NumericVector f(n);
  if(is_ind_null){
    for(int i=0;i<n;++i){
      f[i]=Apply<Rcpp::Matrix<14>::Column,NumericVector,oper,func >(x.column(i),y);
    }
  }else{
    IntegerVector indd(ind);
    for(int i=0;i<n;++i){
      f[i]=Apply<Rcpp::Matrix<14>::Column,NumericVector,oper,func >(x.column(indd[i]-1),y);
    }
  }
  return f;
}

template<Binary_Function oper,Binary_Function func>
double apply_eachrow_helper(SEXP x,SEXP y){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP mat=Rf_duplicate(x);
  double *xx=REAL(mat),*end=xx+ncol*nrow,*yy=REAL(y),y3,*x3,s=0;
  for(;xx!=end;++yy){
    y3=*yy;
    for(x3=xx,xx+=nrow;x3!=xx;++x3){
      s=func(s,oper(*x3,y3));
    }
  }
  return s;
}

template<Binary_Function oper,class T,class RETURN_TYPE,int type>
SEXP eachrow_helper(SEXP x,SEXP y){
  int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
  SEXP mat=PROTECT(Rf_allocMatrix(type,nrow,ncol));
  T *xx=(T *) DATAPTR(x),*xend=xx+ncol*nrow,*yy=(T *) DATAPTR(y),yvalue,*x3;
  RETURN_TYPE *m=(RETURN_TYPE*)DATAPTR(mat);
  for(;xx!=xend;++yy){
    yvalue=*yy;
    for(x3=xx,xx+=nrow;x3!=xx;++x3,++m){
      *m=oper(*x3,yvalue);
    }
  }
  UNPROTECT(1);
  return mat;
}

template<Binary_Function oper,Binary_Function func>
double sum_x_op_y(SEXP x,SEXP y){
  const int n=LENGTH(x);
  double s=0,*startx=REAL(x),*endx=startx+n,*starty=REAL(y);
  for(;startx!=endx;++startx,++starty)
    s=func(s,oper(*startx,*starty));
  return s;
}

template<Binary_Function oper,Binary_Function func>
double sum_x_op_x(SEXP x){
  const int n=LENGTH(x);
  double s=0,*startx=REAL(x),*endx=startx+n;
  for(;startx!=endx;++startx)
    s=func(s,oper(*startx,*startx));
  return s;
}

//sort_unique,len_sort_unique,sort_int
template<typename T>
void max_neg_pos(T *start,T *end,T &mx,T &mn,bool &pos,bool &neg){
  mn=mx=*start;
  T x;
  for(;start!=end;++start){
    x=*start;
    if(x<0){
      neg=true;
      if(mn>x){
        mn=x;
      }
    }else{
      pos=true;
      if(mx<x){
        mx=x;
      }
    }
  }
}

template<Mfunction<int,int,int> Cond,Mfunction<int,int,int> Oper>
int Apply_helper(int *start,int *end,int& val){
    int t=0;
    for(;start!=end;++start){
        if(Cond(*start,val)){
            t=Oper(t,*start);
        }
    }
    return t;
}

template<Binary_Function less_or_greater,Binary_Function min_or_max> // less or greater
NumericVector negative_or_positive(NumericVector &x){
    double v,val=x[0];
    for(auto xx=x.begin()+1;xx!=x.end();++xx){
        v=*xx;
        if(less_or_greater(v,0)){
            if(min_or_max(v,val)){
                val=v;
            }
        }
    }
    return NumericVector::create(val);
}

template<Binary_Function less_or_greater> // less or greater
NumericVector negative_or_positive_min_max(NumericVector &x){
    double v,mn=x[0],mx=mn;
    for(auto xx=x.begin()+1;xx!=x.end();++xx){
        v=*xx;
        if(less_or_greater(v,0)){
            if(v<mn){
                mn=v;
            }else if(v>mx){
                mx=v;
            }
        }
    }
    return NumericVector::create(mn,mx);
}

template<class T>
double nth_simple(T& x,const int& elem){
  nth_element(x.begin(),x.begin()+elem-1,x.end());
  return x[elem-1];
}

template<class T>
double nth_simple(T& x,const int& elem,const bool& descend){
  descend ?
      nth_element(x.begin(),x.begin()+elem-1,x.end(),[&](double a,double b){return a>b;})
  :
      nth_element(x.begin(),x.begin()+elem-1,x.end());

  return x[elem-1];
}

template<class T>
double nth_na_rm(T& x,const int& elem,const bool& descend){
    const int new_end=remove_if(x.begin(),x.end(),R_IsNA)-x.begin();
    if(elem<new_end){
	    descend ?
	    nth_element(x.begin(),x.begin()+elem-1,x.begin()+new_end,[&](double a,double b){return a>b;})
	    :
	    nth_element(x.begin(),x.begin()+elem-1,x.begin()+new_end);
    }
    return x[elem-1];
}

template<class T>
double nth_helper(T& x,const int elem,const bool descend,const bool na_rm){
  return na_rm ? nth_na_rm<T>(x,elem,descend) : nth_simple<T>(x,elem,descend);
}


template<class T>
int nth_index_simple(T& x,const int& elem,const bool& descend){
    IntegerVector ind=seq(1,x.size());
    descend ?
    nth_element(ind.begin(),ind.begin()+elem-1,ind.end(),[&](int i,int j){return x[i-1]>x[j-1];})
        :
        nth_element(ind.begin(),ind.begin()+elem-1,ind.end(),[&](int i,int j){return x[i-1]<x[j-1];});
    
    return ind[elem-1];
}

template<class T>
int nth_index_na_rm(T& x,const int& elem,const bool& descend){
    const int new_end=remove_if(x.begin(),x.end(),R_IsNA)-x.begin();
    IntegerVector ind= seq(1,new_end);
    descend ?
    nth_element(ind.begin(),ind.begin()+((elem<new_end) ? elem-1-new_end : elem-1),ind.end(),[&](int i,int j){return x[i-1]>x[j-1];})
        :
        nth_element(ind.begin(),ind.begin()+((elem<new_end) ? elem-1-new_end : elem-1),ind.end(),[&](int i,int j){return x[i-1]<x[j-1];});
    
    return ind[elem-1];
}

template<class T>
int nth_helper_index(T& x,const int elem,const bool descend,const bool na_rm){
    return na_rm ? nth_index_na_rm<T>(x,elem,descend) : nth_index_simple<T>(x,elem,descend);
}

template<class T>
T nth_simple_n_elems(T& x,const int& elem,const bool& descend){
  descend ?
      nth_element(x.begin(),x.begin()+elem-1,x.end(),[&](double a,double b){return a>b;})
  :
      nth_element(x.begin(),x.begin()+elem-1,x.end());

  return x(span(0,elem-1));
}

template<class T>
T nth_na_rm_n_elems(T& x,const int& elem,const bool& descend){
    const int new_end=remove_if(x.begin(),x.end(),R_IsNA)-x.begin();
    if(elem<new_end){
      descend ?
      nth_element(x.begin(),x.begin()+elem-1,x.begin()+new_end,[&](double a,double b){return a>b;})
      :
      nth_element(x.begin(),x.begin()+elem-1,x.begin()+new_end);
    }
    return x(span(0,elem-1));
}

template<class T>
T nth_helper_n_elems(T& x,const int elem,const bool descend,const bool na_rm){
  return na_rm ? nth_na_rm_n_elems<T>(x,elem,descend) : nth_simple_n_elems<T>(x,elem,descend);
}


template<class T>
T nth_index_simple_n_elems(T& x,const int& elem,const bool& descend){
    vec ind= linspace(1,x.size(),x.size());
    descend ?
    nth_element(ind.begin(),ind.begin()+elem-1,ind.end(),[&](int i,int j){return x[i-1]>x[j-1];})
        :
        nth_element(ind.begin(),ind.begin()+elem-1,ind.end(),[&](int i,int j){return x[i-1]<x[j-1];});
    
    return ind(span(0,elem-1));
}

template<class T>
T nth_index_na_rm_n_elems(T& x,const int& elem,const bool& descend){
    const int new_end=remove_if(x.begin(),x.end(),R_IsNA)-x.begin();
    vec ind= linspace(1,new_end,new_end);
    descend ?
    nth_element(ind.begin(),ind.begin()+((elem<new_end) ? elem-1-new_end : elem-1),ind.end(),[&](int i,int j){return x[i-1]>x[j-1];})
        :
        nth_element(ind.begin(),ind.begin()+((elem<new_end) ? elem-1-new_end : elem-1),ind.end(),[&](int i,int j){return x[i-1]<x[j-1];});
    
    return ind(span(0,elem-1));
}

template<class T>
T nth_helper_index_n_elems(T& x,const int elem,const bool descend,const bool na_rm){
    return na_rm ? nth_index_na_rm_n_elems<T>(x,elem,descend) : nth_index_simple_n_elems<T>(x,elem,descend);
}

#endif
