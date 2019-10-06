\name{Many tests for the dispersion parameter in Poisson distribution}
\alias{colpoisdisp.tests}
\alias{colpois.tests}
\title{
Many tests for the dispersion parameter in Poisson distribution
}

\description{
Many tests for the dispersion parameter in Poisson distribution.
}

\usage{
colpoisdisp.tests(y, alternative = "either", logged = FALSE)
colpois.tests(y, logged = FALSE)
}

\arguments{
\item{y}{
A numerical matrix with count data, 0, 1,...
}
\item{alternative}{
Do you want to test specifically for either over or underspirsion ("either"), overdispersion ("over") or undersispersion ("under")? 
}
\item{logged}{
Set to TRUE if you want the logarithm of the p-value.
}
}

%\details{
%
%}

\value{
A matrix with two columns, the test statistic and the (logged) p-value.
}
\references{
Yang Zhao, James W. Hardin, and Cheryl L. Addy. (2009). A score test for overdispersion in Poisson regression based on the generalized Poisson-2 model. 
Journal of statistical planning and inference 139(4):1514-1521.

Dimitris Karlis and Evdokia Xekalaki (2000). A Simulation Comparison of Several Procedures for 
Testing the Poisson Assumption. Journal of the Royal Statistical Society. Series D (The Statistician), 
49(3): 355-382.

Bohning, D., Dietz, E., Schaub, R., Schlattmann, P. and Lindsay, B. (1994) The distribution of the likelihood 
ratio for mixtures of densities from the one-parameter exponential family. Annals of the Institute of 
Statistical Mathematics, 46(): 373-388.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{poisson.mle}, \link{negbin.mle}, \link{poisson.anova}, \link{poisson.anovas}, \link{poisson_only}
}
}

\examples{
y <- matrix(rnbinom(100* 50, 10, 0.6), ncol = 50)
a1 <- colpoisdisp.tests(y, "over")
b1 <- colpois.tests(y)

y <- matrix(rpois(100* 50, 10), ncol = 50)
a2 <- colpoisdisp.tests(y, "either")
b2 <- colpois.tests(y)
y <- NULL
}

\keyword{ Poisson distribution }
\keyword{ dispersion parameter }
