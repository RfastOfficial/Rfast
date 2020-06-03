\name{Column-wise MLE of some univariate distributions}
\alias{colexpmle}
\alias{colexp2.mle}
\alias{colgammamle}
\alias{colinvgauss.mle}
\alias{collaplace.mle}
\alias{collindley.mle}
\alias{colmaxboltz.mle}
\alias{colnormal.mle}
\alias{colpareto.mle}
\alias{colrayleigh.mle}
\alias{colvm.mle}
\alias{colweibull.mle}
\alias{colnormlog.mle}

\title{
Column-wise MLE of some univariate distributions
}

\description{
Column-wise MLE of some univariate distributions.
}

\usage{
colexpmle(x)
colexp2.mle(x)
colgammamle(x, tol = 1e-07)
colinvgauss.mle(x)
collaplace.mle(x)
collindley.mle(x)
colmaxboltz.mle(x)
colnormal.mle(x)
colpareto.mle(x)
colrayleigh.mle(x)
colvm.mle(x, tol = 1e-07)
colweibull.mle(x, tol = 1e-09, maxiters = 100, parallel = FALSE)
colnormlog.mle(x)
}

\arguments{
\item{x}{
A numerical matrix with data. Each column refers to a different vector of observations of the same distribution.
For exponential, 2 parameter exponential, Weibull, gamma, inverse Gaussian, Maxwell-Boltzman, Lindley, Rayleigh and Pareto
distributions, the numbers must be greater than zero. For the Poisson and geometric distributions, the numbers
must be integers, 0, 1, 2,... For the Normal and Laplace distribution the numbers can take any value. The von
Mises distribution takes values beween 0 and 2 * pi (radians).
}
\item{tol}{
The tolerance value to terminate the Newton-Fisher algorithm.
}
\item{maxiters}{
The maximum number of iterations to implement.
}
\item{parallel}{
Do you want to calculations to take place in parallel? The default value is FALSE
}
}

\details{
For each column, the same distribution is fitted and its parameter and log-likelihood are computed.
}

\value{
A matrix with two, three or five (for the colnormlog.mle) columns. The first one or the first two contain the parameter(s) of the distribution
and the other columns contain the log-likelihood values.
}

\references{
Kalimuthu Krishnamoorthy, Meesook Lee and Wang Xiao (2015). Likelihood ratio tests for comparing several
gamma distributions. Environmetrics, 26(8):571-583.

N.L. Johnson, S. Kotz \& N. Balakrishnan (1994). Continuous Univariate Distributions, Volume 1 (2nd Edition).

N.L. Johnson, S. Kotz \& N. Balakrishnan (1970). Distributions in statistics: continuous univariate distributions,
Volume 2

Sharma V. K., Singh S. K., Singh U. \& Agiwal V. (2015). The inverse Lindley distribution: a stress-strength
reliability model with application to head and neck cancer data. Journal of Industrial and Production
Engineering, 32(3): 162-173.
}

\author{
Michail Tsagris and Stefanos Fafalios

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Stefanos Fafalios
<stefanosfafalios@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vm.mle}, \link{poisson.mle}, \link{normal.mle}, \link{gammamle}
}
}

\examples{
x <- matrix(rnorm(1000 * 50), ncol = 50)
a <- colnormal.mle(x)
b <- collaplace.mle(x)
x <- NULL
}

\keyword{ Column-wise }
\keyword{ continuous distributions }
\keyword{ discrete distributions }
