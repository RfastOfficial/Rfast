//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "Rfast.h"
#include <cmath>
#include "reg_lib.h"


using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
NumericMatrix multinom_regs(NumericVector Y0, NumericMatrix X0, const double tol,
                            const bool logged, const bool parallel, const int maxiters){
  int n = X0.nrow(), D = X0.ncol();

  mat x(X0.begin(), n,D,false);
  mat Y = design_matrix_helper<mat,NumericVector>(Y0);

  rowvec m0 = mean(Y);

  vec poia = indexesOfNum(Y,1);
  int poiasize = poia.size();

  double ini = calc_multinom_ini(Y,conv_to<vec>::from(m0));

  Y.shed_col(0);

  n = Y.n_rows;
  int d = Y.n_cols;

  mat id,b10(2,d,fill::zeros),e0;

  rowvec b0(d);
  m0.shed_col(0);

  b0 = log(m0/(1-m0));
  b10.row(0) = b0;
  e0 = Y.each_row()-m0;

  id = create_id_mat(d);

  double exp20 = exp(20);

  int dx2 = 2*d;
  NumericMatrix ret(D,2);
  colvec one(n,fill::ones);
  /*--------------------------------------------------------*/
    if(parallel){
    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
      mat dera(n,dx2),der2(dx2,dx2,fill::zeros), b1,b2(2,d),m1,m,e1,crossress,X(n,2),xCrossx;
      vec der(dx2),idcoli,idcolj,slv;
      mat::iterator slvit, b2it, b1it;
      int i=0,j=0,ij=0, k=0;
      X.col(0) = one;
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for(int l = 0; l < D; l++) {
        X.col(1) = x.col(l);
        xCrossx = cross_x_y<mat,mat,vec>(X, X);
        for(i = 0; i < d; i++) {
          idcoli = id.col(i);

          dera.col(idcoli[0]) = e0.col(i)%X.col(0);
          dera.col(idcoli[1]) = X.col(1)%e0.col(i);
          for (j = i; j < d; j++) {
            if (i != j) {
              idcolj = id.col(j);
              crossress = -(m0(i) * m0(j))*xCrossx;
              der2(idcolj[0], idcoli[0]) =  crossress(0,0);
              der2(idcolj[0], idcoli[1]) =  crossress(0,1);
              der2(idcolj[1], idcoli[0]) =  crossress(1,0);
              der2(idcolj[1], idcoli[1]) =  crossress(1,1);

              der2(idcoli[0], idcolj[0]) =  crossress(0,0);
              der2(idcoli[0], idcolj[1]) =  crossress(0,1);
              der2(idcoli[1], idcolj[0]) =  crossress(1,0);
              der2(idcoli[1], idcolj[1]) =  crossress(1,1);
            }
            else {
              crossress = ((m0(i) * (1 - m0(i)))) * xCrossx;

              der2(idcoli[0], idcoli[0]) =  crossress(0,0);
              der2(idcoli[0], idcoli[1]) =  crossress(0,1);
              der2(idcoli[1], idcoli[0]) =  crossress(1,0);
              der2(idcoli[1], idcoli[1]) =  crossress(1,1);
            }
          }
        }
        der = conv_to<vec>::from(sum(dera));

        b1 = b10;

        slv = solve(der2, der);

        b2it = b2.begin();
        b1it = b1.begin();
        slvit = slv.begin();
        for(k=0;k<dx2;k++,b1it++,slvit++,b2it++)
          *b2it = (*b1it)+(*slvit);

        ij=2;
        while(ij++<maxiters && sum_with<abs,mat>(b2-b1) > tol) {
          b1 = b2;

          m1 = clamp(exp(X*b1),0,exp20);

          m = m1.each_col()/ (sum(m1,1) + 1);

          e1 = Y - m;
          for(i = 0; i<d; i++) {
            idcoli = id.col(i);
            dera.col(idcoli[0]) = e1.col(i)%X.col(0);
            dera.col(idcoli[1]) = e1.col(i)%X.col(1);

            for (j = 0; j<d; j++) {
              if (i != j) {
                idcolj = id.col(j);

                crossress = -cross_x_y<mat,mat,vec>(X.each_col()%(m.col(i) % m.col(j)), X);
                der2(idcoli[0], idcolj[0]) =  crossress(0,0);
                der2(idcoli[0], idcolj[1]) =  crossress(0,1);
                der2(idcoli[1], idcolj[0]) =  crossress(1,0);
                der2(idcoli[1], idcolj[1]) =  crossress(1,1);

                der2(idcolj[0], idcoli[0]) =  crossress(0,0);
                der2(idcolj[0], idcoli[1]) =  crossress(0,1);
                der2(idcolj[1], idcoli[0]) =  crossress(1,0);
                der2(idcolj[1], idcoli[1]) =  crossress(1,1);
              }
              else {
                crossress = cross_x_y<mat,mat,vec>(X.each_col()%(m.col(i) % (one - m.col(i))), X);

                der2(idcoli[0], idcoli[0]) =  crossress(0,0);
                der2(idcoli[0], idcoli[1]) =  crossress(0,1);
                der2(idcoli[1], idcoli[0]) =  crossress(1,0);
                der2(idcoli[1], idcoli[1]) =  crossress(1,1);
              }
            }
          }
          der = conv_to<vec>::from(sum(dera));

          slv = solve(der2, der);
          b2it = b2.begin();
          b1it = b1.begin();
          slvit = slv.begin();
          for(k=0;k<dx2;k++,b1it++,slvit++,b2it++)
            *b2it = (*b1it)+(*slvit);
        }

        m1.insert_cols(0,one);
        m1 = m1.each_col() / sum(m1,1);

        ret(l,0) = 2 * calcSumLog(m1,poia,poiasize) - ini;
        ret(l,1) = R::pchisq(ret(l,0), d, false, logged);
      }
      #ifdef _OPENMP
      }
      #endif
  }
  else{
    mat dera(n,dx2),der2(dx2,dx2,fill::zeros), b1(2,d),b2(2,d),m1,m,e1,crossress,X(n,2),xCrossx;
    vec der(dx2),idcoli,idcolj,slv;
    mat::iterator slvit, b2it, b1it;
    int i=0,j=0,ij=0, k=0;
    X.col(0) = one;

    for(int l = 0; l < D; l++) {
      X.col(1) = x.col(l);
      xCrossx = cross_x_y<mat,mat,vec>(X, X);
      for(i = 0; i < d; i++) {
        idcoli = id.col(i);

        dera.col(idcoli[0]) = e0.col(i)%X.col(0);
        dera.col(idcoli[1]) = X.col(1)%e0.col(i);
        for (j = i; j < d; j++) {
          if (i != j) {
            idcolj = id.col(j);
            crossress = -(m0(i) * m0(j))*xCrossx;
            der2(idcolj[0], idcoli[0]) =  crossress(0,0);
            der2(idcolj[0], idcoli[1]) =  crossress(0,1);
            der2(idcolj[1], idcoli[0]) =  crossress(1,0);
            der2(idcolj[1], idcoli[1]) =  crossress(1,1);

            der2(idcoli[0], idcolj[0]) =  crossress(0,0);
            der2(idcoli[0], idcolj[1]) =  crossress(0,1);
            der2(idcoli[1], idcolj[0]) =  crossress(1,0);
            der2(idcoli[1], idcolj[1]) =  crossress(1,1);
          }
          else {
            crossress = ((m0(i) * (1 - m0(i)))) * xCrossx;

            der2(idcoli[0], idcoli[0]) =  crossress(0,0);
            der2(idcoli[0], idcoli[1]) =  crossress(0,1);
            der2(idcoli[1], idcoli[0]) =  crossress(1,0);
            der2(idcoli[1], idcoli[1]) =  crossress(1,1);
          }
        }
      }
      der = conv_to<vec>::from(sum(dera));

      b1 = b10;

      slv = solve(der2, der);

      b2it = b2.begin();
      b1it = b1.begin();
      slvit = slv.begin();
      for(k=0;k<dx2;k++,b1it++,slvit++,b2it++)
        *b2it = (*b1it)+(*slvit);

      ij=2;
      while(ij++<maxiters && sum_with<abs,mat>(b2-b1) > tol) {
        b1 = b2;

        m1 = clamp(exp(X*b1),0,exp20);

        m = m1.each_col()/ (sum(m1,1) + 1);

        e1 = Y - m;
        for(i = 0; i<d; i++) {
          idcoli = id.col(i);
          dera.col(idcoli[0]) = e1.col(i)%X.col(0);
          dera.col(idcoli[1]) = e1.col(i)%X.col(1);

          for (j = 0; j<d; j++) {
            if (i != j) {
              idcolj = id.col(j);

              crossress = -cross_x_y<mat,mat,vec>(X.each_col()%(m.col(i) % m.col(j)), X);

              der2(idcoli[0], idcolj[0]) =  crossress(0,0);
              der2(idcoli[0], idcolj[1]) =  crossress(0,1);
              der2(idcoli[1], idcolj[0]) =  crossress(1,0);
              der2(idcoli[1], idcolj[1]) =  crossress(1,1);

              der2(idcolj[0], idcoli[0]) =  crossress(0,0);
              der2(idcolj[0], idcoli[1]) =  crossress(0,1);
              der2(idcolj[1], idcoli[0]) =  crossress(1,0);
              der2(idcolj[1], idcoli[1]) =  crossress(1,1);
            }
            else {
              crossress = cross_x_y<mat,mat,vec>(X.each_col()%(m.col(i) % (one - m.col(i))), X);

              der2(idcoli[0], idcoli[0]) =  crossress(0,0);
              der2(idcoli[0], idcoli[1]) =  crossress(0,1);
              der2(idcoli[1], idcoli[0]) =  crossress(1,0);
              der2(idcoli[1], idcoli[1]) =  crossress(1,1);
            }
          }
        }
        der = conv_to<vec>::from(sum(dera));

        slv = solve(der2, der);
        b2it = b2.begin();
        b1it = b1.begin();
        slvit = slv.begin();
        for(k=0;k<dx2;k++,b1it++,slvit++,b2it++)
          *b2it = (*b1it)+(*slvit);
      }

      m1.insert_cols(0,one);
      m1 = m1.each_col() / sum(m1,1);

      ret(l,0) = 2 * calcSumLog(m1,poia,poiasize) - ini;
      ret(l,1) = R::pchisq(ret(l,0), d, false, logged);
    }
  }

  return ret;
}

RcppExport SEXP Rfast_multinom_regs(SEXP Y0SEXP,SEXP X0SEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP parallelSEXP,SEXP maxitersSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type Y0(Y0SEXP);
    traits::input_parameter< NumericMatrix >::type X0(X0SEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const bool >::type logged(loggedSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    __result = multinom_regs(Y0,X0,tol,logged,parallel,maxiters);
    return __result;
END_RCPP
}
