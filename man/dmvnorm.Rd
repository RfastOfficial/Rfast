\name{Density of the multivariate normal and t distributions}
\alias{dmvnorm}
\alias{dmvt}
\title{
Density of the multivariate normal and t distributions
}
\description{
Density of the multivariate normal and t distributions.
}

\usage{
dmvnorm(x, mu, sigma, logged = FALSE) 
dmvt(x, mu, sigma, nu, logged = FALSE) 
}

\arguments{
\item{x}{
A numerical matrix with the data. The rows correspond to observations and the columns to variables.
}
\item{mu}{
The mean vector.
}
\item{sigma}{
The covariance matrix.
}
\item{nu}{
The degrees of freedom for the multivariate t distribution.
}
\item{logged}{
Should the logarithm of the density be returned (TRUE) or not (FALSE)? 
}
}

\details{
The (log) density of the multivariate normal distribution is calculated for given mean vector and covariance matrix.
}

\value{
A numerical vector with the density values calculated at each vector (row of the matrix x).
}

\references{
Kanti V. Mardia, John T. Kent and John M. Bibby (1979). Multivariate analysis. Academic Press, London.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{rmvnorm}, \link{rmvt}, \link{mvnorm.mle}, \link{iag.mle}
}
}

\examples{
x <- matrnorm(100, 20)
mu <- colmeans(x)
s <- cova(x)
a1 <- dmvnorm(x, mu, s) 
a2 <- dmvt(x, mu, s, 1)
x <- NULL 
}

\keyword{ Multivariate normal distribution }
\keyword{ multivariate t distribution }
\keyword{ density values }
