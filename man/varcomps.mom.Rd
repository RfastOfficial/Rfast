\name{Moment and maximum likelihood estimation of variance components}
\alias{rint.mle}
\alias{varcomps.mom}
\alias{varcomps.mle}
\title{
Moment and maximum likelihood estimation of variance components
}

\description{
Moment and maximum likelihood estimation of variance components.
}

\usage{
rint.mle(x, ina, ranef = FALSE, tol = 1e-09, maxiters = 100)
varcomps.mom(x, ina)
varcomps.mle(x, ina, tol = 1e-09)
}

\arguments{
\item{x}{
A numerical vector with the data.
}
\item{ranef}{
Should the random effects be returned as well? The default value is FALSE.
}
\item{ina}{
A numerical vector with 1s, 2s, 3s and so one indicating the two groups. Be careful, the function is desinged to 
accept numbers greater than zero. Alternatively it can be a factor variable.
}
\item{tol}{
The tolerance level to terminate the golden ratio search. the default value is 10^(-9).
}
\item{maxiters}{
The maximum number of iterations Newton-Raphson will implement.
}
}

\details{
Note that the "varcomps.mle" and "varcomp.mom" work for \bold{balanced designs only}, i.e. for each subject the same number of
measurements have been taken. The "rint.mle" works for both the balanced and unbalanced designs.

The variance components, the variance of the between measurements and the 
variance of the within are estimated using moment estimators. The "colvarcomsp.mom" 
is the moment analogue of a random effects model which uses likelihood estimation ("colvarcomps.mle"). 
It is much faster, but can give negative variance of the random effects, in which case it becomes zero.
 
The maximum likelihood version is a bit slower (try youselves to see the difference), but statistically 
speaking is to be preferred when small samples are available. The reason why it is only a little bit slower 
and not a lot slower as one would imagine is because we are using a closed formula to calculate the two variance
components (Demidenko, 2013, pg. 67-69). Yes, there are closed formulas for linear mixed models.     
}

\value{
For the "varcomps.mom":
A vector with 5 elemets, The MSE, the estimate of the between variance, 
the variance components ratio and a 95\% confidence for the ratio.

For the "varcomps.mle":
a list with a single component called "info". That is a matrix with 3 columns, The MSE, the estimate of 
the between variance and the log-likelihood value. \bold{If ranef = TRUE} a list including "info" 
and an extra component called "ranef" containing the random effects. It is a matrix with the same number of columns as the data. 
Each column contains the randome effects of each variable. 
}

\references{
D.C. Montgomery (2001). Design and analysis of experiments (5th Edition). New York: John Wiley \& Sons. 

Charles S. Davis (2002). Statistical methods for the analysis of repeated measures. New York: Springer-Verlag.

Demidenko E. (2013). Mixed Models: Theory and Applications with R 2nd Edition). New Jersey: John Wiley \& Sons (Excellent book).
}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{colvarcomps.mle}, \link{rint.reg}, \link{rint.regbx}
}
}

\examples{
## example from Montgomery, pages 514-517
x <- c(98,97,99,96,91,90,93,92,96,95,97,95,95,96,99,98)
ina <- rep(1:4, each = 4)
res<-varcomps.mom(x, ina)
res<-varcomps.mle(x, ina)
}

\keyword{ Variance components }
\keyword{ moments estimation }
\keyword{ maximum likelihood estimation }
