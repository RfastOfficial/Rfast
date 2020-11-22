\name{MLE of the tobit model}
\alias{tobit.mle}
\title{
MLE of the tobit model
}

\description{
MLE of the tobit model.
}

\usage{
tobit.mle(y, tol = 1e-09) 
}

\arguments{
\item{y}{
A vector with positive valued data and zero values. If there are no zero values, a simple normal model is fitted in the end. 
}
\item{tol}{
The tolerance level up to which the maximisation stops; set to 1e-09 by default.
}
}

\details{
The tobin model is useful for (univariate) positive data with left censoring at zero. There is the assumption of a latent variable. Tthe values of that variable which are positive
concide with the observed values. If some values are negative, they are left censored and the observed values are zero. 
Instead of maximising the log-likelihood via a numerical optimiser we have used a Newton-Raphson algorithm which is faster. 
}

\value{
A list with three elements including 
\item{iters}{The number of iterations required for the Newton-Raphson to converge.}
\item{loglik}{The value of the maximised log-likelihood.}
\item{param}{The vector of the parameters.}
}

\references{
Tobin James (1958). Estimation of relationships for limited dependent variables. Econometrica. 26(1):24--36. 

https://en.wikipedia.org/wiki/Tobit_model
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{gammamle}, \link{normal.mle}
}
}

\examples{
x <- rnorm(300, 3, 5)
x[ x < 0 ] <- 0   ## left censoring. Values below zero become zero
system.time( for (i in 1:100) tobit.mle(x) )
}

\keyword{ Tobit model }
\keyword{ censored observations }
\keyword{ left censoring }

