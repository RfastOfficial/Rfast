################################
#### Spatial sign covariance matrix
#### Tsagris Michail 2/2015
#### References: A Durre, D Vogel, DE Tyler (2014)
#### The spatial sign covariance matrix with unknown location
#### Journal of Multivariate Analysis, 130: 107--117.
#### http://arxiv.org/pdf/1307.5706v2.pdf
#### mtsagris@yahoo.gr
################################
#[export]
sscov <- function(x, me = NULL, tol = 1e-09) {
  n <- dim(x)[1]  ## sample size
  if ( is.null(me) )  me <- Rfast::spat.med(x, tol)  ## spatial median of x
  y <- Rfast::eachrow(x, me, oper = "-")
  rs <- sqrt ( Rfast::rowsums(y^2) )
  crossprod( y / rs ) / n  ## SSCM
}
