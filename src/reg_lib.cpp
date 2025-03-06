//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "reg_lib.h"
#include "Rfast/templates.h"
#include "my_k_sorted_array.h"
#include <algorithm>
#include "mn.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

//[[Rcpp::plugins(openmp)]]

double average_value(vec y, a_node* my_ar, const int size){
  double sum=0.0;
  a_node* it = my_ar;
  for(int i=0;i<size;i++,it++){
    sum+= y(it->index);
  }

  return sum/size;
}

double weighted_average_value(vec y, a_node* my_ar, const int size){
  double sum=0.0, divider = 0.0,acos;
  a_node* it = my_ar;
  for(int i=0;i<size;i++,it++){
    acos = std::exp(std::acos(-it->value));
    sum+= y(it->index)/acos;
    divider+=1/acos;
  }

  return sum/divider;
}

double weighted_most_frequent_value(vec y, a_node* my_ar, const int size){
  std::map<int,double> counts;
  a_node* tmp = my_ar;
  double acos;

  for(int i=0;i<size;i++,tmp++){
      acos = std::exp(std::acos(-tmp->value));
      counts[(int)y(tmp->index)]+=1/acos;
  }

  map<int, double>::iterator tmpit;
  int mostFrequent = -1;
  double maxCount = 0;
  for(tmpit = counts.begin(); tmpit!= counts.end(); tmpit++) {
    if(tmpit->second > maxCount) {
      mostFrequent = tmpit->first;
      maxCount = tmpit->second;
    }
  }
  return mostFrequent;
}

double most_frequent_value(vec y, a_node* my_ar, const int size){
  std::map<int,int> counts;
  a_node* tmp = my_ar;
  for(int i=0;i<size;i++,tmp++){
    counts[(int)y(tmp->index)]++;
  }

  map<int, int>::iterator tmpit;
  int mostFrequent = -1;
  int maxCount = 0;
  for(tmpit = counts.begin(); tmpit!= counts.end(); tmpit++) {
    if(tmpit->second > maxCount) {
      mostFrequent = tmpit->first;
      maxCount = tmpit->second;
    }
  }
  return mostFrequent;
}

double vmf_mle2(double nR, const int n, const double tol, const double maxiters){
  double apk, k2, R = nR/n, R2 = R*R, k1 = R * (2 - R2)/(1 - R2);

  int i = 2;
  if(k1 < 1e+05){
    apk = R::bessel_i(k1,1, 1)/R::bessel_i(k1,0, 1);
    k2 = k1 - (apk - R)/(1 - apk*apk - 1/k1 * apk);
    while (i++<maxiters && abs(k2 - k1) > tol) {
      k1 = k2;
      apk = R::bessel_i(k1,1,1)/R::bessel_i(k1,0,1);
      k2 = k1 - (apk - R)/(1 - apk*apk - 1/k1 * apk);
    }
    return k2;
  }

  return k1;
}
//lik1 <-  - 0.5 * sum( mu^2 ) + sum( log1p( tau * ptau * con / exp(f * tau^2) ) )
double calc_spml_loglik(mat::col_iterator mu1, mat::col_iterator mu2, vec::iterator tau, vec::iterator ptau, const int size){
  double f= -0.5, con = 2.506628274631;
  double ret1 = 0.0, ret2 = 0.0;

  for(int i = 0;i<size;i++,ptau++,tau++,mu2++,mu1++){
    ret1+= (*mu1)*(*mu1)+(*mu2)*(*mu2);
    ret2+=  log1p(((*tau)*(*ptau))*con/ exp(f*(*tau)*(*tau)));
  }

  return -0.5*ret1+ret2;
}

double spml_mle2(mat u, vec ci2, vec cisi, vec si2, const int n, const double tol, const int maxiters){
  vec su(2);
  su(0) = sum(u.col(0)), su(1) = sum(u.col(1));

  double nR = sqrt(sum(su%su)), kappa = vmf_mle2(nR, n, tol, maxiters);

  vec mu = su*(1/nR);

  vec mu1 = mu*kappa;
  double f = -0.5, con = 2.506628274631;
  vec tau = u*mu1, ptau = pnormc(tau);

  vec rat = ptau/(exp(f * tau%tau)/con + tau % ptau);

  vec psit = tau + rat;
  vec psit2 = 2 - rat%(tau + rat);

  vec der(2);
  der[0] = sum(u.col(0)%psit) - n * mu1[0];
  der[1] = sum(u.col(1)%psit) - n * mu1[1];

  double dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;

  double down = dera2 * derb2 - derab*derab;
  vec mu2(2);
  mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
  mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

  int i = 2;
  while (i++<maxiters && sum(abs(mu2 - mu1)) > tol) {
    mu1 = mu2;
    tau = u*mu1;
    ptau = pnormc(tau);
    rat = ptau/(exp(f * (tau%tau))/con + tau % ptau);
    psit = tau + rat;
    psit2 = 2 - rat%(tau + rat);
    der[0] = sum(u.col(0)%psit) - n * mu1[0];
    der[1] = sum(u.col(1)%psit) - n * mu1[1];
    dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;
    down = dera2 * derb2 - derab*derab;
    mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
    mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

  }

  return -0.5 * n * (mu2[0]*mu2[0]+mu2[1]*mu2[1]) + sum_with<log1p, colvec>(((tau % ptau) * con)/exp(f*tau%tau)) - n * 1.83787706640935;
}

void my_pow2(vec inp,double *out,const double power,const int sz){
  for(double *startx=&inp[0],*starty=out,*end=startx+sz;startx!=end;++startx,++starty){
    *starty=std::pow(*startx,power);
  }

  return;
}

double calc_f(vec nix, double n, vec ni2hi2, double S, double x, int size){
  double sum1 = 0.0, sum2 = 0.0;

  for(int i = 0; i < size; i++){
    sum1+=log1p(nix[i]);
    sum2+=ni2hi2[i]/(1+nix[i]);
  }

  return sum1+n*log(S-x*sum2);
}

vec gold_rat3(double n, vec ni, vec ni2, double S, vec hi2,const int size, const double tol=1e-07){
  double a = 0, b = 50;
  const double ratio=0.618033988749895;
  double x1=b-ratio*b, x2=ratio*b;
  vec nix1 = ni*x1, nix2 = ni*x2, ni2hi2 = ni2%hi2;

  double f1 = calc_f(nix1, n, ni2hi2, S, x1, size);
  double f2 = calc_f(nix2, n, ni2hi2, S, x2, size);
  double bmina = b - a;
  while (abs(bmina)>tol){
    if(f2>f1){
      b=x2;
      bmina = b - a;
      x2=x1;
      f2=f1;
      x1=b - ratio * (bmina);
      nix1 = ni*x1;
      f1 = calc_f(nix1, n, ni2hi2, S, x1, size);
    }
    else {
      a=x1;
      bmina = b - a;
      x1=x2;
      f1=f2;
      x2=a + ratio * (bmina);
      nix2 = ni*x2;
      f2 = calc_f(nix2, n, ni2hi2, S, x2, size);
    }
  }
  vec ret(2);
  ret(0) = 0.5*(x1+x2);
  ret(1) = (f1+f2)/2;

  return ret;
}

colvec log1pColvec(colvec input, int sz){
  colvec ret(sz);
  colvec::iterator iter1 = input.begin(),iter2 = ret.begin(), end = input.end();

  for(; iter1 != end; ++iter1,++iter2){
    *iter2 = log1p(*iter1);
  }

  return ret;
}

//lik = n * anew[0] + anew[1] * tmpsumX - sum(y1% log1pColvec(exp( anew[0] + anew[1] * tmpX),n));
double get_geom_lik(const double a1,const double a2, const double sX, double *tmpX, double *y1,const int n){
  double ret = 0.0;
  for(int i =0;i<n;i++,y1++,tmpX++){
    ret+= (*y1)*std::log1p(std::exp(a1+a2*(*tmpX)));
  }
  return n*a1+a2*sX-ret;
}

vec indexesOfNum(mat m, int num){
  int sz = m.n_cols*m.n_rows;
  vec tmp(sz);
  int i,j = 0;

  for(i=0; i<sz;i++){
    if(m(i)==num){
      tmp(j++)=i;
    }
  }

    tmp.resize(j);

    return tmp;
}

mat create_id_mat(int d){
  mat ret(2,d);
  ret(0,0) = 0;
  ret(1,0) = 1;

  for(int i=1;i<d;i++){
    ret(0,i) = ret(0,i-1)+2;
    ret(1,i) = ret(1,i-1)+2;
  }
  return ret;
}

double calc_multinom_ini(mat Y1,vec m0){
  double ret = 0.0;
  int n=Y1.n_rows;
  rowvec logm0 = conv_to<rowvec>::from(log(m0));
  for(int i = 0;i<n;i++){
    ret+=sum(Y1.row(i)%logm0);
  }

  return 2*ret;
}

double calcSumLog(mat ma, vec poia, int sz){
  double ret = 0.0;
  for(int i=0; i < sz;i++){
    ret+=log(ma(poia[i]));
  }
  return ret;
}

mat colvec_mat_cbind(vec v, mat m){
  int n = m.n_rows,d=m.n_cols;
  mat ret(n,d+1);
  ret.col(0)= v;

  for(int i =1; i<d+1;i++){
    ret.col(i) = m.col(i-1);
  }

  return ret;
}

double calcDevRes(mat p,vec y,mat est){
  int psize = p.n_rows;
  double summ =0.0;

  for(int i=0;i<psize;i++){
    if(y(i)==1){
      if(p(i,0) == 0){
        summ += est(i,0);
      }
      else{
        summ+=log(p(i,0));
      }
    }
    else{
      if(p(i,0) == 1){
        summ+= est(i,0);
      }
      else{
        summ+=log(1-p(i,0));
      }
    }
  }

  return summ;
}

vec subvec(vec v, int start, int size){
  vec ret(size);
  for(int i=start,j = 0;i<size+start;i++,j++){
    ret(j) = v(i);
  }
  return ret;
}

double varcomps_mle2(vec x, IntegerVector ina,int n,const double tol) {
  //Author: Manos Papadakis
  const int N=x.size(),d=N/n;
  vec y=abs(x-mean(x)),syina=group_sum_helper<vec,vec,IntegerVector>(y,ina,nullptr,&n);
  double sy2=sum(sqrt(syina)), a=0,ratio=2.0/(sqrt(5) + 1),sy=sum(sqrt(y)),b=sy/N,s=b;
  double x1=b-ratio*b,x2=ratio*b;
  double se=s-x1;
  double f1=N*log(se)+n*log1p(d*x1/se)+sy/se-x1/(se*se+d*x1*se)*sy2;
  se=s-x2;
  double f2=N*log(se)+n*log1p(d*x2/se)+sy/se-x2/(se*se+d*x2*se)*sy2;

  while (abs(b-a)>tol){
    if(f2>f1){
      b=x2;
      x2=x1;
      f2=f1;
      x1=b - ratio * (b - a);
      se=s - x1;
      f1=N * log(se) + n * log1p(d * x1 / se) + sy/se - x1 / (se*se + d * x1 * se) * sy2 ;
    } else {
      a=x1;
      x1=x2;
      f1=f2;
      x2=a + ratio * (b - a);
      se=s - x2;
      f2=N * log(se) + n * log1p(d * x2 / se) + sy/se - x2 / (se*se + d * x2 * se) * sy2;
    }
  }

  return -0.5 * f2 - N*0.5 * 1.837877;
}

mat varcomps_mle3(vec x, IntegerVector ina,const int N, int n,const bool ranef, const double tol,const int maxiters) {
  const int d=N/n;
  int i = 2;
  vec y=abs(x-mean(x)),syina=group_sum_helper<vec,vec,IntegerVector>(y,ina,nullptr,&n);

  double sy2=sum(sqrt(syina)), a=0,ratio=2.0/(sqrt(5) + 1),sy=sum(sqrt(y)),b=sy/N,s=b;
  double x1=b-ratio*b,x2=ratio*b;
  double se=s-x1;
  double f1=N*log(se)+n*log1p(d*x1/se)+sy/se-x1/(se*se+d*x1*se)*sy2;
  se=s-x2;
  double f2=N*log(se)+n*log1p(d*x2/se)+sy/se-x2/(se*se+d*x2*se)*sy2;

  while (i++ < maxiters && abs(b-a)>tol){
    if(f2>f1){
      b=x2;
      x2=x1;
      f2=f1;
      x1=b - ratio * (b - a);
      se=s - x1;
      f1=N * log(se) + n * log1p(d * x1 / se) + sy/se - x1 / (se*se + d * x1 * se) * sy2 ;
    } else {
      a=x1;
      x1=x2;
      f1=f2;
      x2=a + ratio * (b - a);
      se=s - x2;
      f2=N * log(se) + n * log1p(d * x2 / se) + sy/se - x2 / (se*se + d * x2 * se) * sy2;
    }
  }

  mat ret;
  if(ranef){
    if(n>4){
      ret = mat(2,n);
    }
    else{
      ret = mat(2,4);
    }
  }
  else{
    ret = mat(1,4);
  }

  const double tau=(a+b)/2.0;

  ret(0,0) = i-1;
  ret(0,2) = tau;
  ret(0,1) = s - tau;
  ret(0,3) = -0.5 * (f2 + N * 1.837877);


  if(ranef){
    ret.row(1) = conv_to<rowvec>::from(ret(0,2)/(ret(0,2) + ret(0,1)/d) * syina/d);
  }

  return ret;
}

vec weibull_mle2(vec x, int n, const double tol, const int maxiters){
  int i=2;
  vec lx = log(x),lx2 = lx%lx, y = x;
  double mlx = sum(lx)/n, co = sum(y%lx),sy = sum(y), fb = 1+mlx-co/sy,fb2 = -1 -(sum(y%lx2)*sy-co*co)/(sy*sy);
  double b1 = 1,b2 = 1 - fb/fb2;

  while (++i<maxiters && sum(abs(b2 - b1)) > tol) {
    b1 = b2;
    my_pow2(x,&y[0],b1,n);
    co = sum(y % lx);
    sy = sum(y);
    fb = 1/b1 + mlx - co/sy;
    fb2 = -1/(b1*b1) - (sum(y % lx2) * sy - co*co)/(sy*sy);
    b2 = b1 - fb/fb2;
  }

  vec param(3);
  double theta = pow(sy/n,1/b2);
  my_pow2(conv_to<vec>::from(x/theta),&y[0],b2,n);
  param[0] = b2;
  param[1] = theta;
  param[2] = n * log(b2) - n * b2 * log(theta) + (b2 - 1) * n * mlx - sum(y);

  return param;
}
