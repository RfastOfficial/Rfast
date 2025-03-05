#define ARMA_64BIT_WORD
#include <vector>
#include <string>
#include <algorithm>
#include "mn.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

int proper_size(int nrw, int ncl)
{
  return ncl * (ncl - 1) * 0.5;
}

double sum_max_elems(colvec x, colvec y)
{
  double maxs = 0.0;
  for (unsigned int i = 0; i < x.n_elem; ++i)
  {
    maxs += max(x[i], y[i]);
  }
  return maxs;
}

double sum_min_elems(colvec x, colvec y)
{
  double mins = 0.0;
  for (unsigned int i = 0; i < x.n_elem; ++i)
  {
    mins += min(x[i], y[i]);
  }
  return mins;
}

colvec max_elems(colvec x, colvec y)
{
  colvec maxs(x.n_elem,fill::none);
  for (unsigned int i = 0; i < x.n_elem; ++i)
  {
    maxs[i] = max(x[i], y[i]);
  }
  return maxs;
}

mat colMaxElems(mat x, colvec y)
{
  mat maxs(x.n_rows,x.n_cols,fill::none);
  for (unsigned int i = 0; i < x.n_cols; ++i)
  {
    maxs.col(i) = max_elems(x.col(i),y);
  }
  return maxs;
}

rowvec euclidean_norm(mat &x)
{
	return sqrt(sum(square(x), 0));
}

bool my_compare_order_second(const pr<double,int>& a,const pr<double,int>& b){
  return a.second<b.second;
}

rowvec operator/(colvec x,double s){
  rowvec f(x.n_elem);
  for(rowvec::iterator ff=f.begin(),xx=x.begin();ff!=f.end();++ff,++xx){
    *ff=*xx/s;
  }
  return f;
}

//regression
double regression_only_col(colvec x, colvec &y) {
  int n=x.size();
  double SSO=var(y)*(double)(n-1),SS1=0.0,F1=0.0;
  mat z(n,2,fill::ones),tr_z(2,n);
  colvec b(2);
  vec res(n);
  z.col(1)=x;
  tr_z=z.t();
  b=inv(tr_z*z)*tr_z*y;
  res=y-z*b;
  SS1=var(res)*(n-1);
  F1=(SSO/SS1-1)*(n-2);
  return F1;
}

//diri_nr_type2,gamma
double trigamma ( double x)
{
  double a = 0.0001;
  double b = 5.0;
  double b2 =  0.1666666667;
  double b4 = -0.03333333333;
  double b6 =  0.02380952381; 
  double b8 = -0.03333333333;
  double value;
  double y;
  double z;
  
  z = x;
  //
  //  Use small value approximation if X <= A.
  //
  if ( x <= a )
  {
    value = 1.0 / x / x;
    return value;
  }
  //
  //  Increase argument to ( X + I ) >= B.
  //
  value = 0.0;
  
  while ( z < b )
  {
    value = value + 1.0 / z / z;
    z = z + 1.0;
  }
  //
  //  Apply asymptotic formula if argument is B or greater.
  //
  y = 1.0 / z / z;
  
  value = value + 0.5 *
    y + ( 1.0
            + y * ( b2
                      + y * ( b4
                      + y * ( b6
                      + y *   b8 )))) / z;
                      
                      return value;
}

//diri_nr_type2,gamma
double digamma(double x) {
  double result = 0, xx, xx2, xx4;
  for ( ; x < 7; ++x){
    result -= 1/x;
  }
  x -= 1.0/2.0;
  xx = 1.0/x;
  xx2 = xx*xx;
  xx4 = xx2*xx2;
  result += log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
  return result;
}

//floyd
void i4mat_floyd( int n, NumericVector &a ){
  int i,j,k;
  const double i4_huge = 2147483647;
  for ( k = 0; k < n; k++ ){
    for ( j = 0; j < n; j++ ){
      if ( a[k+j*n] < i4_huge ){
        for ( i = 0; i < n; i++ ){
          if ( a[i+k*n] < i4_huge ){
            a[i+j*n] = std::min( a[i+j*n], a[i+k*n] + a[k+j*n] );
          }
        }
      }
    }
  }
}

void i4mat_floyd_with_paths( const int n, NumericVector &a,NumericVector &p ){
  int i,j,k;
  const double i4_huge = 2147483647;
  for ( k = 0; k < n; k++ ){
    for ( j = 0; j < n; j++ ){
      if ( a[k+j*n] < i4_huge ){
        for ( i = 0; i < n; i++ ){
          if ( a[i+k*n] < i4_huge ){
            a[i+j*n] = std::min( a[i+j*n], a[i+k*n] + a[k+j*n] );
            p[i+j*n] = k;
          }
        }
      }
    }
  }
}


//comb_n
void combn(arma::vec& vals, const int n, const unsigned int start_idx, 
    std::vector<double>& combn_data, double*& combn_col) {
  if (!n) {
    for (unsigned int i = 0; i < combn_data.size(); ++i) {
      *combn_col++ = combn_data[i];
    }
    return;
  }
  for (unsigned int i = start_idx; i <= (vals.size() - n); ++i) {
    combn_data.at(combn_data.size() - n) = vals(i);
    combn(vals, n - 1, i + 1, combn_data, combn_col);
  }
}

//rmdp
uvec Order_rmdp(colvec& x){
  uvec ind=linspace<uvec>(0,x.n_elem-1,x.n_elem);
  std::stable_sort(ind.begin(),ind.end(),[&](int i,int j){return x[i]<x[j];});
  return ind;
}

//rmdp
rowvec colvar_rmdp(mat& x){
  rowvec nyr1=x.row(0),nyr2=x.row(1);
  return 0.5*(arma::square(nyr1) + arma::square(nyr2)) - nyr1%nyr2;
}

//dists
double sum_pow(colvec x,const double p){
  const int sz=x.size();
  double s=0;
  for(double *startx=&x[0],*end=startx+sz;startx!=end;++startx)
    s+=std::pow(*startx,p);
  return s;
}

//Design_matrix
umat design_matrix_helper_big(CharacterVector x) {
  int i=0;
  const int n=x.size();
  CharacterVector tmp=sort_unique(x);
  CharacterVector::iterator xx=x.begin(),leksi_bg,leksi_en;
  umat Final(n,tmp.size(),fill::zeros);
  for(leksi_bg=tmp.begin(),leksi_en=tmp.end(),i=0;xx!=x.end();++xx,++i)
    Final(i,std::lower_bound(leksi_bg,leksi_en,*xx)-leksi_bg)=1;
  return Final;
}

//varcomps_mle
NumericVector minus_mean(NumericVector& x,const double k){
  NumericVector y(x.size());
  for(NumericVector::iterator xx=x.begin(),yy=y.begin();x.end()-xx;++xx,++yy){
    *yy=*xx-k;
  }
  return y;
}

//vecdist
void minus_c(double f[],double &x,double *y,int offset,int &len){
  double *ff=f;
  for(int i=0;i<len;++i,ff+=offset,++y){
    *ff=std::abs(x-*y);
  }
}

//squareform,Round
int my_round(const double x){
  return ((int(x)*10)%10>4) ? int(x)+1 : x;
}

static long double powers[] = {0,1e+1,1e+2,1e+3,1e+4,1e+5,1e+6,1e+7
								,1e+8,1e+9,1e+10,1e+11,1e+12,1e+13
								,1e+14,1e+15,1e+16};

//Round
double my_round_gen_na_rm(double x,const int& dg){
  if(R_IsNA(x)){
    return x;
  }
  long long int t=powers[dg+1];
  const bool nx=x<0;
  long long int y= nx ? -x*t : x*t;
  const int m=y%10;
  y= m>4 ? y+10-m : y-m ;
  x=y;
  return nx ? -x/t : x/t ; 
}


double my_round_gen_simple(double x,const int& dg){
  long long int t=powers[dg+1];
  const bool nx=x<0;
  long long int y= nx ? -x*t : x*t;
  const int m=y%10;
  y= m>4 ? y+10-m : y-m ;
  x=y;
  return nx ? -x/t : x/t ; 
}

//Norm
double sumsqr(NumericMatrix &x){
  double s=0,v;
  for(double *start=x.begin(),*end=x.end();start!=end;++start){
    v=*start;
    s+=v*v;
  }
  return std::sqrt(s);    
}

//col/row True
int True(int *start,int *end){
  int t=0;
  for(;start!=end;++start){
    if(*start){
      ++t;
    }
  }
  return t;
}

//all
bool my_all(int* start,int *end){
  for(;start!=end;++start){
    if(!(*start)){
      return false;
    }
  }
  return true;
}

//any
bool my_any(int* start,int *end){
  for(;start!=end;++start){
    if(*start){
      return true;
    }
  }
  return false;
}

//spml_mle
colvec pnormc(colvec x){
  for(double *xx=&x[0],*endx=&x[x.n_elem];xx!=endx;++xx){
    *xx=R::pnorm(*xx,0,1,1,0);
  }
  return x;
}

//spml_mle
double sum_abs(mat x,mat y){
  double s=0;
  for(unsigned int i=0;i<x.n_elem;++i){
    s+=std::abs(x[i]-y[i]);
  }
  return s;
}

//hash2lists
NumericVector toNumbers(string x,const string spliter){
  NumericVector f;
  x+=spliter;
  const char *split=spliter.c_str();
  char *xx=(char*)x.c_str();
  char *token = std::strtok(xx, split);
  while (token != nullptr) {
    f.push_back(std::atof(token));
    token = std::strtok(nullptr, split);
  }
  return f;
}

//bincomb
IntegerVector combine(IntegerVector x,IntegerVector y){
  const int n=x.size(),p=y.size(),z=n+p;
  IntegerVector f(z);
  f[Range(0,n-1)]=x;
  f[Range(n,z-1)]=y;
  return f;
}

//dista
icolvec get_k_indices(rowvec x,const int& k){
  icolvec ind=linspace<icolvec>(1,x.size(),x.size());
  std::sort(ind.begin(),ind.end(),[&](icolvec::elem_type i,icolvec::elem_type j){return x[i-1]<x[j-1];});
  return ind(span(0,k-1));
}

colvec get_k_values(rowvec x, const int &k)
{
	sort(x.begin(), x.end());
	return conv_to<colvec>::from(x.subvec(0, k - 1));
}

bool check_if_is_finite(double x)
{
	return x > 0 and !R_IsNA(x);
}

double calcDevRes(colvec p,colvec y,colvec expyhat){
  int psize = p.n_elem;
  double s=0.0;
  for(int i=0;i<psize;i++){
    if(y(i)==1){
      if(p(i) == 0){
        s+= expyhat(i);
      }
      else{
        s+=log(p(i));
      }
    }
    else{
      if(p(i) == 1){
        s+= expyhat(i);
      }
      else{
        s+=log(1-p(i));
      }
    }
  }
  
  return s;
}
