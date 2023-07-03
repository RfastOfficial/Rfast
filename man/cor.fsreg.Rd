\name{Correlation based forward regression}
\alias{cor.fsreg}
\title{
Correlation based forward regression.
}

\description{
Correlation based forward regression.
}

\usage{
cor.fsreg(y, x, ystand = TRUE, xstand = TRUE, threshold = 0.05, 
tolb = 2, tolr = 0.02, stopping = "BIC") 
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A matrix with data, the predictor variables. 
}
\item{ystand}{
If this is TRUE the response variable is centered. The mean is subtracted from every value.
}
\item{xstand}{
If this is TRUE the independent variables are standardised. 
}
\item{threshold}{
The significance level, set to 0.05 by default. Bear in mind that the logarithm of it is used, as the 
logarithm of the p-values is calculated at every point. This will avoid numerical overflows and
small p-values, less than the machine epsilon, being returned as zero.
}
\item{tolb}{
If we see only the significane of the variables, many may enter the linear regression model. For this reason, 
we also use the BIC as a way to validate the inclusion of a candidate variable. If the BIC difference between 
two successive models is less than the tolerance value, the variable will not enter the model, even if it 
statistically significant. Set it to 0 if you do not want this extra check. 
}
\item{tolr}{
This is an alternative to the BIC change and it uses the adjusted coefficient of determination. If the increase in the 
adjusted \eqn{R^2} is more than the tolr continue.
}
\item{stopping}{
This refers to the type of extra checking to do. If you want the BIC check, set it to "BIC". If you want the adjusted \eqn{R^2} check 
set this to "ar2". Or, if you want both of them to take place, both of these criteria to be satisfied make this "BICR2".
}
}

\details{
The forward regression tries one by one the variables using the F-test, basically partial F-test every time for 
the latest variable. This is the same as testing the significance of the coefficient of this latest enetered 
variable. Alternatively the correlation can be used and this case the partial correlation coefficient. There is 
a direct relationship between the t-test statistic and the partial correlation coefficient. Now, instead of 
having to calculate the test statistic, we calculate the partial correlation coefficient. Using Fisher's 
z-transform we get the variance imediately. The partial correlation coefficient, using Fisher's z-transform, 
and the partial F-test (or the coefficient's t-test statistic) are not identical. They will be identical for 
large sample sizes though.
}

\value{
A matrix with three columns, the index of the selected variables, the logged p-value and the the test 
statistic value and the BIC or adjusted \eqn{R^2} of each model. In the case of stopping="BICR2" both of 
these criteria will be returned.
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
\code{ \link{score.glms}, \link{univglms}, \link{logistic_only},  \link{poisson_only}, \link{regression}
}
}

\examples{
## 200 variables, hence 200 univariate regressions are to be fitted
x <- matrnorm(200,  100)
y <- rnorm(200)
cor.fsreg(y, x)
x <- NULL
}

\keyword{ Forward regression }
\keyword{ Pearson correlation }
\keyword{ partial correlation }