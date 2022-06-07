\name{Quasi Poisson regression for count data}
\alias{qpois.reg}
\alias{qpois.regs}
\title{
Quasi Poisson regression
}

\description{
Quasi Poisson regression.
}

\usage{
qpois.reg(x, y, full = FALSE, tol = 1e-09,maxiters = 100) 
qpois.regs(x, y, tol = 1e-09, logged = FALSE) 
}

\arguments{
\item{x}{
For the "qpois.reg" a matrix with data, the predictor variables. This can be a matrix or a 
data frame. For the "qpois.regs" this must be a numerical matrix, where each columns denotes 
a variable. 
}
\item{y}{
A numerical vector with positive discrete data. 
}
\item{full}{
If this is FALSE, the coefficients, the deviance and the estimated phi parameter will be 
returned only. If this is TRUE, more information is returned. 
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. This is set to \eqn{10^{-9}} by default.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)? 
}
\item{maxiters}{
The maximum number of iterations before the Newton-Raphson is terminated automatically.
}
}

\details{
We are using the Newton-Raphson, but unlike R's built-in function "glm" we do no checks 
and no extra calculations, or whatever. Simply the model, unless the user requests for the 
Wald tests of the coefficients.  The "qpois.regs" is to be used for very many univariate 
regressions. The "x" is a matrix in this case and the significance of each variable 
(column of the matrix) is tested. 
}

\value{
For the "prop.reg" a list including:
When full is FALSE 
\item{be}{
The regression coefficients.
}
\item{devi}{
The deviance of the model.
}
\item{varb}{
The covariance matrix of the beta coefficients.
}
\item{phi}{
The phi parameter, the estimate of dispersion.
}
When full is TRUE, the additional item is:
\item{info}{
The regression coefficients, their standard error, their Wald test statistic and their p-value. 
}
For the "prop.regs" a two-column matrix with the test statistics (Wald statistic) and 
the associated p-values (or their loggarithm). 
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989.
}

\author{
Manos Papadakis and Marios Dimitriadis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> and 
Marios Dimitriadis <kmdimitriadis@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{prop.reg} \link{univglms}, \link{score.glms}, \link{poisson_only}
}
}

\examples{
\donttest{
y <- rnbinom(100, 10, 0.6)
x <- matrix(rnorm(100*3), ncol = 3)
mod1 <- glm(y ~ x, quasipoisson)
res<-summary(mod1)
res<-qpois.reg(x, y, full = TRUE)
res<-qpois.regs(x, y)
}
}

\keyword{ GLMs }
\keyword{ Quasi Poisson regression }
