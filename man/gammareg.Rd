\name{Gamma regression with a log-link}
\alias{gammareg}
\alias{gammacon}
\title{
Gamma regression with a log-link
}

\description{
Gamma regression with a log-link.
}

\usage{
gammareg(y, x, tol = 1e-07, maxiters = 100)
gammacon(y, tol = 1e-08, maxiters =50)
}

\arguments{
\item{y}{
The dependent variable, a numerical variable with non negative numbers.
}
\item{x}{
A matrix or data.frame with the indendent variables.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations that can take place in the regression.
}
}

\details{
The gamma.reg fits a Gamma regression with a log-link. The gamma.con fits a Gamma regression
with a log link with the intercept only ( glm(y ~ 1, Gamma(log) ) ). 
}

\value{
A list including:
\item{deviance}{
The deviance value.
}
\item{phi}{
The dispersion parameter (\eqn{\phi}) of the regression. This is necessary if you want to 
perform an F hypothesis test for the significance of one or more independent variables.
}
\item{be}{
The regression coefficient(s). 
}
\item{info}{
The number of iterations, the deviance and the dispersion parameter.
}
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{ \link{gammaregs}, \link{normlog.reg}, \link{invgauss.reg} }
}

\examples{
y <- abs( rnorm(100) )
x <- matrix( rnorm(100 * 2), ncol = 2)
mod <- glm(y ~ x, family = Gamma(log) )
res<-summary(mod)
\donttest{
res<-gammareg(y, x)
}
mod <- glm(y ~ 1, family = Gamma(log) )
res<-summary(mod)
res<-gammacon(y)
}

