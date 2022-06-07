\name{Circular or angular regression}
\alias{spml.reg}
\title{
Circular or angular regression
}

\description{
Regression with circular dependent variable and Euclidean or categorical independent variables.
}

\usage{
spml.reg(y, x, tol = 1e-07, seb = FALSE, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, it can be a numerical vector with data expressed in radians or it can be a matrix with
two columns, the cosinus and the sinus of the circular data. The benefit of the matrix is that if the function 
is to be called multiple times with the same response, there is no need to transform the vector every time into 
a matrix. 
}
\item{x}{
The independent variable(s). Can be Euclidean or categorical (factor variables).
}
\item{tol}{
The tolerance value to terminatate the Newton-Raphson algorithm.
}
\item{seb}{
Do you want the standard error of the estimates to be returned? TRUE or FALSE.
}
\item{maxiters}{
The maximum number of iterations to implement.
}
}

\details{
The Newton-Raphson algorithm is fitted in this regression as described in Presnell et al. (1998).
}

\value{
A list including:
\item{iters}{
The number of iterations required until convergence of the EM algorithm.
}
\item{be}{
The regression coefficients.
}
\item{seb}{
The standard errors of the coefficients.
}
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{seb}{
The covariance matrix of the beta values.
}
}

\references{
Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for 
directional data. Journal of the American Statistical Association, 93(443): 1068-1077.
}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis
<papadakm95@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spml.mle}, \link{iag.mle}, \link{acg.mle}
}
}

\examples{
\donttest{
x <- rnorm(100)
z <- cbind(3 + 2 * x, 1 -3 * x)
y <- cbind( rnorm(100,z[ ,1], 1), rnorm(100, z[ ,2], 1) )
y <- y / sqrt( rowsums(y^2) )
a1 <- spml.reg(y, x)
y <- atan( y[, 2] / y[, 1] ) + pi * I(y[, 1] < 0) 
a2 <- spml.reg(y, x)
}
}

\keyword{ Circular regression }
\keyword{ projected normal }
