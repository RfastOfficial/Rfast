\name{MLE of (hyper-)spherical distributions}
\alias{vmf.mle}
\alias{multivmf.mle}
\alias{acg.mle}
\alias{iag.mle}
\title{
MLE of (hyper-)spherical distributions
}

\description{
MLE of (hyper-)spherical distributions.
}

\usage{
vmf.mle(x, tol = 1e-07)
multivmf.mle(x, ina, tol = 1e-07, ell = FALSE)
acg.mle(x, tol = 1e-07)
iag.mle(x, tol = 1e-07)
}

\arguments{
\item{x}{
A matrix with directional data, i.e. unit vectors.
}
\item{ina}{
A numerical vector with discrete numbers starting from 1, i.e. 1, 2, 3, 4,... or a factor variable. Each number denotes a sample or group. 
If you supply a continuous valued vector the function will obviously provide wrong results.
}
\item{ell}{
This is for the multivmf.mle only. Do you want the log-likelihood returned? The default value is TRUE. 
}
\item{tol}{
The tolerance value at which to terminate the iterations.  
}
}

\details{
For the von Mises-Fisher, the normalised mean is the mean direction. For the concentration parameter, a Newton-Raphson is implemented. 
For the angular central Gaussian distribution there is a constraint on the estimated covariance matrix; its trace is equal to the number of variables. 
An iterative algorithm takes place and convergence is guaranteed. Newton-Raphson for the projected normal distribution, on the sphere, is implemented as well. Finally, 
the von Mises-Fisher distribution for groups of data is also implemented. 
}

\value{
For the von Mises-Fisher a list including:
\item{loglik}{
The maximum log-likelihood value.
}
\item{mu}{
The mean direction.
}
\item{kappa}{
The concentration parameter.
}
For the multi von Mises-Fisher a list including:
\item{loglik}{
A vector with the maximum log-likelihood values if ell is set to TRUE. Otherwise NULL is returned.
}
\item{mi}{
A matrix with the group mean directions.
}
\item{ki}{
A vector with the group concentration parameters.
}
For the angular central Gaussian a list including:
\item{iter}{
The number if iterations required by the algorithm to converge to the solution.
}
\item{cova}{
The estimated covariance matrix.
}
For the spherical projected normal a list including:
\item{iters}{
The number of iteration required by the Newton-Raphson. 
}
\item{mesi}{
A matrix with two rows. The first row is the mean direction and the second is the mean vector. 
The first comes from the second by normalising to have unit length.
}
\item{param}{
A vector with the elements, the norm of mean vector, the log-likelihood and the log-likelihood of the spherical uniform distribution. 
The third value helps in case you want to do a log-likleihood ratio test for uniformity.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra, S. (2012). A short note on parameter approximation for von Mises-Fisher distributions: and a fast implementation of Is(x). 
Computational Statistics, 27(1): 177--190.

Tyler D. E. (1987). Statistical analysis for the angular central Gaussian distribution on the sphere.
Biometrika 74(3): 579-589.

Paine P.J., Preston S.P., Tsagris M and Wood A.T.A. (2017). An Elliptically Symmetric Angular Gaussian Distribution. Statistics and Computing (To appear). 
}

\author{
Michail Tsagris
R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{racg}, \link{vm.mle}, \link{rvmf}
}
}

\examples{
m <- c(0, 0, 0, 0)
s <- cov(iris[, 1:4])
x <- racg(100, s)
mod <- acg.mle(x)
mod
res<-cov2cor(mod$cova)  ## estimated covariance matrix turned into a correlation matrix
res<-cov2cor(s)  ## true covariance matrix turned into a correlation matrix
res<-vmf.mle(x)
x <- rbind( rvmf(100,rnorm(4), 10), rvmf(100,rnorm(4), 20) )
a <- multivmf.mle(x, rep(1:2, each = 100) )
}

\keyword{ von Mises-Fisher distribution }
\keyword{ central angular Gaussian distribution }
\keyword{ maximum likelihood estimation }
\keyword{ projected normal distribution }
\keyword{ directional data }
\keyword{ spherical data }