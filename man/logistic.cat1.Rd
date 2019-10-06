\name{Logistic or Poisson regression with a single categorical predictor}
\alias{logistic.cat1}
\alias{poisson.cat1}
\title{
Logistic or Poisson regression with a single categorical predictor
}

\description{
Logistic or Poisson regression with a single categorical predictor.
}

\usage{
logistic.cat1(y, x, logged = FALSE)
poisson.cat1(y, x, logged = FALSE) 
}

\arguments{
\item{y}{
A numerical vector with values 0 or 1. 
}
\item{x}{
A numerical vector with discrete numbers or a factor variable. This is suppose to be a categorical predictor. 
If you supply a continuous valued vector the function will obviously provide wrong results. 
\bold{Note:} For the "binomial.anova" if this is a numerical vector it must contain strictly positive numbers, i.e.
1, 2, 3, 4, ..., no zeros are allowed. 
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
There is a closed form solution for the logistic regression in the case of a single predictor variable. 
See the references for more information. 
}

\value{
\item{info}{
A matrix similar to the one produced by the glm command. The estimates, their standard error, the Wald value and
the relevant p-value.
}
\item{devs}{
For the logistic regression case a vector with the null and the residual deviances, their difference and the 
significance of this difference.  
}
\item{res}{
For the Poisson regression case a vector with the log likelihood ratio test statistic value and its significance. 
}
}
\references{
Stan Lipovetsky (2015). Analytical closed-form solution for binary logit regression by categorical
predictors. Journal of Applied Statistics, 42(1): 37--49.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{poisson.anova}, \link{poisson.anovas}, \link{anova}, \link{logistic_only}, \link{poisson_only}
}
}

\examples{
y <- rbinom(20000, 1, 0.6)
x <- as.factor( rbinom(20000, 3, 0.5) )
system.time( a1 <- logistic.cat1(y, x) )
system.time( a2 <- glm(y ~ x, binomial) ) 
a1  ;  a2

y <- rpois(20000, 10)
x <- as.factor( rbinom(20000, 3, 0.5) )
system.time( a1 <- poisson.cat1(y, x) )
system.time( a2 <- glm(y ~ x, poisson) ) 
a1   ;  a2

x<-y<-a1<-a2<-NULL
}

\keyword{ Logistic regression }
\keyword{ Poisson regression }
\keyword{ single categorical predictor }
\keyword{ analysis of variance }