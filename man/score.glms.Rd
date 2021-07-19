\name{Many score based regressions}
\alias{score.glms}
\alias{score.multinomregs}
\alias{score.negbinregs}
\alias{score.weibregs}
\alias{score.betaregs}
\alias{score.gammaregs}
\alias{score.expregs}
\alias{score.invgaussregs}
\alias{score.ztpregs}
\alias{score.geomregs}
\title{
Many score based regressions
}

\description{
Many score based GLM regressions.
}

\usage{
score.glms(y, x, oiko = NULL, logged = FALSE) 
score.multinomregs(y, x, logged = FALSE) 
score.negbinregs(y, x, logged = FALSE) 
score.weibregs(y, x, logged = FALSE) 
score.betaregs(y, x, logged = FALSE) 
score.gammaregs(y, x, logged = FALSE) 
score.expregs(y, x, logged = FALSE)
score.invgaussregs(y, x, logged = FALSE)
score.ztpregs(y, x, logged = FALSE)
score.geomregs(y, x, logged = FALSE)
}

\arguments{
\item{y}{
A vector with either discrete or binary data for the Poisson or negative binomial and binary logistic regression respectively. 
A vector with discrete values or factor values for the multinomial regression. If the vector 
is binary and choose multinomial regression the function checks and transfers to the binary logistic regression. 

For the Weibull, gamma and exponential regressions they must be strictly positive data, 
lifetimes or durations for example. For the beta regression they must be numbers between 0 and 1. For the 
zero truncated Poisson regression (score.ztpregs) they must be integer valued data strictly greater than 0. 
}
\item{x}{
A matrix with data, the predictor variables. 
}
\item{oiko}{
This can be either "poisson" or "binomial". If you are not sure leave it NULL and the function will check 
internally. 
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
}

\details{
Instead of maximising the log-likelihood via the Newton-Raphson algorithm in order to perform the hypothesis 
testing that \eqn{\beta_i=0} we use the score test. This is dramatcially faster as no model needs to be fitted. 
The first derivative (score) of the log-likelihood is known and in closed form and under the null hypothesis 
the fitted values are all equal to the mean of the response variable y. The variance of the score is also known 
in closed form. The test is not the same as the likelihood ratio test. It is size correct nonetheless but it 
is a bit less efficient and less powerful. For big sample sizes though (5000 or more) the results are the same. 
We have seen via simulation studies is that it is size correct to large sample sizes, at elast a few thousands. 
You can try for yourselves and see that even with 500 the results are pretty close. The score test is pretty faster 
then the classical log-likelihood ratio test.  
}

\value{
A matrix with two columns, the test statistic and its associated p-value. For the Poisson and logistic 
regression the p-value is derived via the t distribution, whereas for the multinomial regressions 
via the \eqn{\chi^2} distribution. 
}

\references{
Tsagris M., Alenazi A. and Fafalios S. (2020). Computationally efficient univariate filtering for massive data. 
Electronic Journal of Applied Statistical Analysis, 13(2):390-412.

Hosmer DW. JR, Lemeshow S. and Sturdivant R.X. (2013). Applied Logistic Regression. New Jersey ,Wiley, 3rd Edition.

Campbell M.J. (2001). Statistics at Square Two: Understand Modern Statistical Applications in Medicine, pg. 112.
London, BMJ Books. 
 
Draper N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.

McCullagh Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.

Agresti Alan (1996). An introduction to categorical data analysis. New York: Wiley.

Joseph M.H. (2011). Negative Binomial Regression. Cambridge University Press, 2nd edition.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{univglms}, \link{logistic_only},  \link{poisson_only}, \link{regression}
}
}

\examples{
x <- matrnorm(500, 500)
y <- rbinom(500, 1, 0.6)   ## binary logistic regression
a2 <- score.glms(y, x)
y <- rweibull(500, 2, 3)   
a <- score.weibregs(y, x) 
mean(a[, 2] < 0.05) 
x <- NULL
}

\keyword{ GLMS }
\keyword{ multinomial regressions }
\keyword{ score based tests }