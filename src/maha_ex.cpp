// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec mahaInt(arma::mat & X,arma::vec & mu,arma::mat & sigma,const bool isChol);

/* 
 *  Internal C++ function for Mahalanobis distance
*/
RcppExport SEXP Rfast_mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP isChol)
{
  using namespace Rcpp;
  
  try{
    
    arma::mat X_ = as<arma::mat>(X);
    arma::vec mu_ = as<arma::vec>(mu);  
    arma::mat sigma_ = as<arma::mat>(sigma); 
    bool isChol_ = as<bool>(isChol);
    
    NumericVector dist = wrap( mahaInt(X_, mu_, sigma_, isChol_) );
    dist.attr( "dim" ) = R_NilValue;
    
    return dist;
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return wrap(NA_REAL);
}