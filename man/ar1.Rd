\name{Estimation of an AR(1) model}
\alias{ar1}
\alias{colar1}
\title{
Estimation of an AR(1) model
}

\description{
Estimation of an AR(1) model.
}

\usage{
ar1(y, method = "cmle") 
colar1(y, method = "cmle")
}

\arguments{
\item{y}{
For the case of \bold{ar1} this is a vector of time series. For the case of \bold{colar1} this is a matrix where weach column represents a time series. 
}
\item{method}{
This can be either "cmle" for conditional maximum likelihood or "yw" for the Yule-Walker equations.
}
}

\details{
Instead of the classical MLE for the AR(1) model which requires numerical optimsation (Newton-Raphson for example) we estimate the parameters of the AR(1) model using 
conditional maximum likelihood. This procedure is described in Chapter 17 in Lee (2006). In some, it assumes that the first observation is deterministic and hence 
conditioning on that observation, there is a closed form solution for the parameters. The second alternative is to use the method of moments and hence the Yule-Walker equations. 
}

\value{
\item{param}{
For the case of \bold{ar1} this is a vector with three elements, the constant term, the \eqn{\phi} term (lag coefficient) and the variance.
For the case of \bold{colar1} this is a matrix with three columns, eahc of which carries the same aforementioned elements. 
}
}

\references{
http://econ.nsysu.edu.tw/ezfiles/124/1124/img/Chapter17_MaximumLikelihoodEstimation.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{rm.lines}, \link{varcomps.mle}, \link{rm.anovas}
}
}

\examples{
y <- as.vector(lh)
ar1(y)
ar(y, FALSE, 1, "ols")

ar1(y, method = "yw")
ar(y, FALSE, 1, "yw")

a1 <- colar1(cbind(y, y) )
b1 <- colar1(cbind(y, y), method = "yw")
}

\keyword{ Time series }
\keyword{ AR(1) model }
\keyword{ conditional MLE }