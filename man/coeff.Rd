\name{Coefficient matrices}
\alias{coeff}
\title{
Coefficient matrices
}

\description{
Coefficient matrices.
}

\usage{
coeff(x, method,vector = FALSE)
}

\arguments{
\item{x}{
A matrix with data. The distances will be calculated between pairs of rows. 
In the case of \bold{vecdist} this is a vector. 
For the \bold{haversine} distance it must be a matrix with two columns, 
the first column is the latitude and the second the longitude. 
}
\item{method}{
See details for the available methods.
}
\item{vector}{
For return a vector instead a matrix.
}
}

\details{
\itemize{
\item bhattacharyya : \eqn{  \sum \sqrt(P_i * Q_i)}
}
}

\value{
A square matrix with the pairwise coefficients. 
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{dista}, \link{Dist}
}
}

\examples{
x <- matrix(rnorm(50 * 10), ncol = 10)
a1 <- coeff(x,"bhattacharyya")
x<-a1<-NULL
}

\keyword{ Coefficient matrices }