\name{Fitting a Dirichlet distribution via Newton-Rapshon}
\alias{diri.nr2}
\title{
Fitting a Dirichlet distribution via Newton-Rapshon
}

\description{
Fitting a Dirichlet distribution via Newton-Rapshon.
}

\usage{
diri.nr2(x, type = 1, tol = 1e-07)
}

\arguments{
\item{x}{
A matrix containing the compositional data. Zeros are not allowed.
}
\item{type}{
Type 1 uses a vectorised version of the Newton-Raphson (Minka, 2012). In high dimensions this is to be preferred. 
If the data are too concentrated, regardless of the dimensions, this is also to be preferrred. 
Type 2 uses the regular Newton-Raphson, with matrix multiplications. In small dimensions this can be considerably
faster. 
}
\item{tol}{
The tolerance level idicating no further increase in the log-likelihood.
}
}

\details{
Maximum likelihood estimation of the parameters of a Dirichlet distribution is performed via Newton-Raphson. 
Initial values suggested by Minka (2012) are used. 
}

\value{
A list including:
\item{loglik}{
The value of the log-likelihood.
}
\item{param}{
The estimated parameters.
}
}

\references{
Minka Thomas (2012). Estimating a Dirichlet distribution. Technical report.

Ng Kai Wang, Guo-Liang Tian, and Man-Lai Tang (2011). Dirichlet and related distributions: Theory, methods and 
applications. John Wiley & Sons.
}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>
}

\seealso{
\code{\link{beta.mle}
}
}

\examples{
x <- matrix( rgamma(100 * 4, c(5, 6, 7, 8), 1), ncol = 4)
x <- x / rowsums(x)
res<-diri.nr2(x) 
}

\keyword{ Dirichlet distribution }
\keyword{ maximum likelihood estimation }
\keyword{ Newton-Raphson }
