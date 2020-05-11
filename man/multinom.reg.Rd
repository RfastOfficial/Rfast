\name{Multinomial regression}
\alias{multinom.reg}
\title{
Multinomial regression
}

\description{
Multinomial regression.
}

\usage{
multinom.reg(y, x, tol = 1e-07, maxiters = 50)
}

\arguments{
\item{y}{
The response variable. A numerical or a factor type vector.
}
\item{x}{
A matrix or a data.frame with the predictor variables.
}
\item{tol}{
This tolerance value to terminate the Newton-Raphson algorithm. 
}
\item{maxiters}{
The maximum number of iterations Newton-Raphson will perform.
}
}

\value{
A list including:
\item{iters}{
The number of iterations required by the Newton-Raphson.
}
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{be}{
A matrix with the estimated regression coefficients.
}
}

\references{
Bohning, D. (1992). Multinomial logistic regression algorithm. Annals of the 
Institute of Statistical Mathematics, 44(1): 197-200.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos 
Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{glm_logistic}, \link{score.multinomregs} \link{logistic_only}
}
}

\examples{
\dontrun{
y <- iris[, 5]
x <- matrnorm(150, 3)
res <- multinom.reg(y, x)
}
}

\keyword{ Multinomial distribution }
\keyword{ regression }