\name{Exponential empirical likelihood for a one sample mean vector hypothesis testing}
\alias{mv.eeltest1}
\title{
Exponential empirical likelihood for a one sample mean vector hypothesis testing
}

\description{
Exponential empirical likelihood for a one sample mean vector hypothesis testing.
}

\usage{
mv.eeltest1(x, mu, tol = 1e-06)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{mu}{
The hypothesized mean vector.
}
\item{tol}{
The tolerance value used to stop the Newton-Raphson algorithm.
}
}

\details{
Multivariate hypothesis test for a one sample mean vector. This is a non parametric test and it works for univariate and multivariate data. The p-value is currently computed only asymptotically (no bootstrap calibration at the moment).
}

\value{
A list including:
\item{p}{
The estimated probabiities.
}
\item{lambda}{
The value of the Lagrangian parameter \eqn{\lambda}.
}
\item{iters}{
The number of iterations required by the newton-Raphson algorithm.
}
\item{info}{
The value of the log-likelihood ratio test statistic along with its corresponding p-value.
}
}

\references{
Jing Bing-Yi and Andrew TA Wood (1996). Exponential empirical likelihood is not Bartlett correctable. Annals of Statistics 24(1): 365-369.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{james}, \link{mv.eeltest2}
}
}

\examples{
x <- Rfast::rmvnorm(100, numeric(10), diag( rexp(10, 0.5) ) )
res<-mv.eeltest1(x, numeric(10) )
}

\keyword{ Hypothesis test }
\keyword{ mean vector }
