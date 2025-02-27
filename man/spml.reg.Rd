\name{Circular or angular regression(s)}
\alias{spml.reg}
\alias{spml.regs}
\title{
Circular or angular regression(s)
}

\description{
Regression with circular dependent variable and Euclidean or categorical independent variables.
}

\usage{
spml.reg(y, x, tol = 1e-07, seb = FALSE, maxiters = 100)
spml.regs(y, x, tol = 1e-07, logged = FALSE, maxiters = 100, parallel = FALSE)
}

\arguments{
\item{y}{
The dependent variable, it can be a numerical vector with data expressed in radians or it can be a matrix with
two columns, the cosinus and the sinus of the circular data. The benefit of the matrix is that if the function 
is to be called multiple times with the same response, there is no need to transform the vector every time into 
a matrix. 
}
\item{x}{
The independent variable(s) of the single regression model. Can be Euclidean or categorical (factor variables).
For the many regression models, a matrix with continuous variables.
}
\item{tol}{
The tolerance value to terminatate the Newton-Raphson algorithm.
}
\item{seb}{
Do you want the standard error of the estimates to be returned? TRUE or FALSE.
}
\item{logged}{
Do you want the logarithm of the p-value be returned? TRUE or FALSE.
}
\item{maxiters}{
The maximum number of iterations to implement.
}
\item{parallel}{
Do you want the regression models to run in parallel? Default value is FALSE.
}
}

\details{
The Newton-Raphson algorithm is fitted in this regression as described in Presnell et al. (1998).

For the case of spml.regs(), for each colum of x a circual regression model is fitted and the 
hypothesis testing of no association between y and this variable is performed.
}

\value{
For the spml.reg() a list including:
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

For the spml.regs() a matrix with two columns, the test statistics and their associated (log) p-values.
}

\references{
Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for 
directional data. Journal of the American Statistical Association, 93(443): 1068-1077.
}

\author{
Michail Tsagris and Stefanos Fafalios.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Stefanos Fafalios
<stefanosfafalios@gmail.com>.
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

x <- matrnorm(100, 5)
a3 <- spml.regs(y, x)
}
}

\keyword{ Circular regression }
\keyword{ projected normal }
