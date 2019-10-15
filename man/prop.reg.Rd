\name{Quasi binomial regression for proportions}
\alias{prop.reg}
\alias{prop.regs}
\title{
Quasi binomial regression for proportions
}

\description{
Quasi binomial regression for proportions.
}

\usage{
prop.reg(y, x, varb = "quasi", tol = 1e-09, maxiters = 100) 
prop.regs(y, x, varb = "quasi", tol = 1e-09, logged = FALSE, maxiters = 100) 
}

\arguments{
\item{y}{
A numerical vector proportions. 0s and 1s are allowed. 
}
\item{x}{
For the "prop.reg" a matrix with data, the predictor variables. This can be a matrix or a data frame.
For the "prop.regs" this must be a numerical matrix, where each columns denotes a variable.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. This is set to \eqn{10^{-9}} by default.
}
\item{varb}{
The type of estimate to be used in order to estimate the covariance matrix of the regression coefficients. 
There are two options, either "quasi" (default value) or "glm". See the references for more information. 
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?  
}
\item{maxiters}{
The maximum number of iterations before the Newton-Raphson is terminated automatically.
}
}

\details{
We are using the Newton-Raphson, but unlike R's built-in function "glm" we do no checks and no extra calculations, 
or whatever. Simply the model. The "prop.regs" is to be used for very many univariate regressions. The "x" is a 
matrix in this case and the significance of each variable (column of the matrix) is tested. The function accepts binary 
responses as well (0 or 1).
}

\value{
For the "prop.reg" function a list including:
\item{iters}{
The number of iterations required by the Newton-Raphson.
}
\item{varb}{
The covariance matrix of the regression coefficients.
}
\item{phi}{
The phi parameter is returned if the input argument "varb" was set to "glm", othwerise this is NULL.
}
\item{info}{
A table similar to the one produced by "glm" with the estimated regression coefficients, their standard error,
Wald test statistic and p-values.
}
For the "prop.regs" a two-column matrix with the test statistics (Wald statistic) and the associated p-values 
(or their loggarithm). 
}

\references{
Papke L. E. & Wooldridge J. (1996). Econometric methods for fractional response variables with 
an application to 401(K) plan participation rates. Journal of Applied Econometrics, 11(6): 619--632.

McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{anova_propreg} \link{univglms}, \link{score.glms}, \link{logistic_only}
}
}

\examples{
y <- rbeta(100, 1, 4)
x <- matrix(rnorm(100 * 3), ncol = 3)
a <- prop.reg(y, x)
y <- rbeta(100, 1, 4)
x <- matrix(rnorm(400 * 100), ncol = 400)
b <- prop.regs(y, x)
res<-mean(b[, 2] < 0.05)
}

\keyword{ GLMs }
\keyword{ Quasi regression }
\keyword{ fractional response }