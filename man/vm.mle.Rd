\name{MLE of some circular distributions}
\alias{vm.mle}
\alias{spml.mle}
\alias{wrapcauchy.mle}
\title{
MLE of some circular distributions
}

\description{
MLE of some circular distributions.
}

\usage{
vm.mle(x, tol = 1e-09)
spml.mle(x, tol = 1e-09, maxiters = 100)
wrapcauchy.mle(x, tol = 1e-09)
}

\arguments{
\item{x}{
A numerical vector with the circular data. They must be expressed in radians. For the "spml.mle" this can also
be a matrix with two columns, the cosinus and the sinus of the circular data.
}
\item{tol}{
The tolerance level to stop the iterative process of finding the MLEs.
}
\item{maxiters}{
The maximum number of iterations to implement.
}
}

\details{
The parameters of the von Mises, the bivariate angular Gaussian and wrapped Cauchy distributions are estimated. For the Wrapped Cauchy,
the iterative procedure described by Kent and Tyler (1988) is used. As for the von Mises distribution,
we use a Newton-Raphson to estimate the concentration parameter. The angular Gaussian is described, in the regression setting in Presnell
et al. (1998).
}

\value{
A list including:
\item{iters}{
The iterations required until convergence. This is returned in the wrapped Cauchy distribution only.
}
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{param}{
A vector consisting of the estimates of the two parameters, the mean direction for both distributions
and the concentration parameter kappa and the rho for the von Mises and wrapped Cauchy respectively.
}
\item{gamma}{
The norm of the mean vector of the angualr Gaussian distribution.
}
\item{mu}{
The mean vector of the angular Gaussian distribution.
}
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra S. (2012). A short note on parameter approximation for von Mises-Fisher distributions:
and a fast implementation of Is(x). Computational Statistics, 27(1): 177-190.

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data.
Journal of the American Statistical Association, 93(443): 1068-1077.

Kent J. and Tyler D. (1988). Maximum likelihood estimation for the wrapped Cauchy distribution.
Journal of Applied Statistics, 15(2): 247--254.
}

\author{
Michail Tsagris and Stefanos Fafalios

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Stefanos Fafalios
<stefanosfafalios@gmail.com>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmf.mle}, \link{rvonmises}, \link{rvmf}
}
}

\examples{
y <- rcauchy(100, 3, 1)
x <- y %% (2 * pi)
res<-vm.mle(x)
res<-spml.mle(x)
res<-wrapcauchy.mle(x)
x <- NULL
}

\keyword{ von Mises distribution }
\keyword{ bivariate angular Gaussian }
\keyword{ wrapped Cauchy distribution }
\keyword{ circular data }
