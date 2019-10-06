\name{Tests for the dispersion parameter in Poisson distribution}
\alias{poisdisp.test}
\alias{pois.test}
\title{
Tests for the dispersion parameter in Poisson distribution
}

\description{
Tests for the dispersion parameter in Poisson distribution.
}

\usage{
poisdisp.test(y, alternative = "either", logged = FALSE)
pois.test(y, logged = FALSE)
}

\arguments{
\item{y}{
A numerical vector with count data, 0, 1,...
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
A vector with two elements, the test statistic and the (logged) p-value. 
}
\references{
Yang Zhao, James W. Hardin, and Cheryl L. Addy. (2009). A score test for overdispersion in Poisson regression 
based on the generalized Poisson-2 model. Journal of statistical planning and inference 139(4): 1514-1521.

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
y <- rnbinom(500, 10, 0.6)
res<-poisdisp.test(y, "either")
res<-poisdisp.test(y, "over")
res<-pois.test(y)

y <- rpois(500, 10)
res<-poisdisp.test(y, "either")
res<-poisdisp.test(y, "over")
res<-pois.test(y)
}

\keyword{ Poisson distribution }
\keyword{ dispersion parameter }
