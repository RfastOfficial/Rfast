\name{BIC (using partial correlation) forward regression}
\alias{bic.corfsreg}
\title{
BIC (using partial correlation) forward regression
}

\description{
BIC (using partial correlation) forward regression.
}

\usage{
bic.corfsreg(y, x, tol = 2) 
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A matrix with data, the predictor variables. 
}
\item{tol}{
If the BIC difference between two successive models is less than the tolerance value, the variable will not enter the model.
}
}

\details{
The forward regression tries one by one the variables using the F-test, basically partial F-test every time for 
the latest variable. This is the same as testing the significance of the coefficient of this latest enetered 
variable. Alternatively the correlation can be used and this case the partial correlation coefficient. There is 
a direct relationship between the t-test statistic and the partial correlation coefficient. Now, instead of 
having to calculate the test statistic, we calculate the partial correlation coefficient. The largest partial correlation
indicates the candidate variable to enter the model. If the BIC of the regression model with that variable included, reduces, 
less than "tol" from the previous model without this variable, the variable enters. 
}

\value{
A matrix with two columns, the index of the selected variable(s) and the BIC of each model. The first line is always 
0 and the BIC of the model with no predictor variables. 
}

\references{
Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cor.fsreg}, \link{score.glms}, \link{univglms}, \link{logistic_only},  
\link{poisson_only}, \link{regression}
}
}

\examples{
## 200 variables, hence 200 univariate regressions are to be fitted
x <- matrix( rnorm(200 * 200), ncol = 200 )
y <- rnorm(200)
system.time( a1 <- bic.corfsreg(y, x) )
system.time( a2 <- cor.fsreg(y, x) )
x <- NULL
}

\keyword{ Forward regression }
\keyword{ BIC }
\keyword{ partial correlation }