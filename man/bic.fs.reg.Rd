\name{BIC forward regression with generalised linear models}
\alias{bic.fs.reg}
\title{
BIC forward regression with generalised linear models
}

\description{
BIC forward regression with generalised linear models.
}

\usage{
bic.fs.reg(y, x, tol = 2, type = "logistic") 
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A matrix with data, the predictor variables. 
}
\item{tol}{
If the BIC difference between two successive models is less than the tolerance value, 
the variable will not enter the model.
}
 	
\item{type}{
If you have a binary dependent variable, put "logistic". If you have count data, put "poisson".
}
}

\details{
The forward regression tries one by one the variables using the BIC at each step for 
the latest variable. If the BIC of the regression model with that variable included, 
is less than "tol" from the previous model without this variable, the variable enters. 
}

\value{
A matrix with two columns, the index of the selected variable(s) and the BIC of each model. 
}

\references{
Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.
}

\author{
Marios Dimitriadis

R implementation and documentation: Marios Dimitriadis <kmdimitriadis@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{fs.reg}, \link{bic.corfsreg}, \link{cor.fsreg}, \link{score.glms}, \link{univglms}, \link{logistic_only},  
\link{poisson_only}, \link{regression}
}
}

\examples{
\dontrun{
x <- matrix(rnorm(200 * 50), ncol = 50)
## 200 variables, hence 200 univariate regressions are to be fitted
y <- rbinom(200, 1, 0.5)
a <- bic.fs.reg(y, x)
x <- NULL
}
}

\keyword{ Forward regression }
\keyword{ BIC }
