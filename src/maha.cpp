/* 
Copyright (C) 2014 Matteo Fasiolo  matteo.fasiolo@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */

/*
 *  Fast computation of Mahalanobis distance
 *
 * See ?maha() for a description of the arguments and output.
 */ 
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/* 
*  Internal C++ function for Mahalanobis distance
*/
Rcpp::NumericVector mahaInt(arma::mat & X,  
                  arma::vec & mu,  
                  arma::mat & sigma,
                  const bool isChol = false)
{
  using namespace arma;
  
  // Some sanity checks 
  if(mu.n_elem != sigma.n_cols) Rcpp::stop("The mean vector has a different dimensions from the covariance matrix.");
  if(X.n_cols != sigma.n_cols)  Rcpp::stop("The number of columns of X is different from the dimension of the covariance matrix.");
  
  // Calculate transposed cholesky dec. unless sigma is alread a cholesky dec.
  mat cholDec;
  if( isChol == false ) {
    cholDec = trimatl(chol(sigma).t());
  }
  else{
    cholDec = trimatl(sigma.t()); 
    if(any(cholDec.diag() <= 0.0))  Rcpp::stop("The supplied cholesky decomposition has values <= 0.0 on the main diagonal.");
  }
  
  vec D = cholDec.diag();

  Rcpp::NumericVector Out(X.n_rows);
  vec out(Out.begin(), Out.size(), false);
  
  // Declaring some private variables
  uint32_t d = X.n_cols;
  uint32_t n = X.n_rows;
  
  vec tmp(d);  
  
  double acc;
  uint32_t icol, irow, ii;  
  
  // For each of the "n" random vectors, forwardsolve the corresponding linear system.
  // Forwardsolve because I'm using the lower triangle Cholesky.
  for(icol = 0; icol < n; icol++)
  {
    
    for(irow = 0; irow < d; irow++)
    {
      acc = 0.0;
      
      for(ii = 0; ii < irow; ii++) acc += tmp.at(ii) * cholDec.at(irow, ii);
      
      tmp.at(irow) = ( X.at(icol, irow) - mu.at(irow) - acc ) / D.at(irow);
    }
    
    out.at(icol) = sum(square(tmp)); 
  }
  
return Out;
}




/* 
  #Equivalent R function:
  
  .fastMahalanobis <- function(X, mean, mcov)
  {
    dec <- chol(mcov)
    tmp <- forwardsolve(t(dec), t(X) - mean )
    colSums( tmp ^ 2 )
  }

*/