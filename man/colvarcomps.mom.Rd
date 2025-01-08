\name{Many moment and maximum likelihood estimations of variance components}
\alias{colvarcomps.mom}
\alias{colvarcomps.mle}
\title{
Many moment and maximum likelihood estimations of variance components
}

\description{
Many moment and maximum likelihood estimations of variance components.
}

\usage{
colvarcomps.mom(x, id, parallel = FALSE) 
colvarcomps.mle(x, id, ranef = FALSE, tol= 1e-08, maxiters = 100, 
parallel = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where each column refers to a different sample of subjects. 
}
\item{id}{
A numerical vector indicating the subject. You must put consecutive numbers and no zero values. 
Alternatively this can be a factor variable.
}
\item{ranef}{
Do you also want the random effects to be returned? TRUE or FALSE.
}
\item{tol}{
The tolerance level to terminate the golden ratio search. 
}
\item{maxiters}{The maximum number of iterations to perform.
}
\item{parallel}{
Should the computations run in parallel? TRUE or FALSE.
}
}

\details{
Note that the "colvarcomp.mom" works for \bold{balanced designs only}, i.e. for each subject the same number of
measurements have been taken. The "colvarcomps.mle" works for unbalanced as well.

The variance components, the variance of the between measurements and the 
variance of the within are estimated using moment estimators. The "colvarcomps.mom" 
is the moment analogue of a random effects model which uses likelihood estimation ("colvarcomps.mle"). 
It is much faster, but can give negative variance of the random effects, in which case it becomes zero.
 
The maximum likelihood version is a bit slower (try youselves to see the difference), but statistically 
speaking is to be preferred when small samples are available. The reason why it is only a little bit slower 
and not a lot slower as one would imagine is because we are using a closed formula to calculate the two variance
components (Demidenko, 2013, pg. 67-69). Yes, there are closed formulas for linear mixed models.     
}

\value{
For the "colvarcomps.mom":
A matrix with 5 columns, The MSE, the estimate of the between variance, 
the variance components ratio and a 95\% confidence 
for the ratio.

For the "colvarcomps.mle":
\bold{If ranef = FALSE} 
a list with a single component called "info". That is a matrix with 3 columns, The MSE, the estimate of 
the between variance and the log-likelihood value. \bold{If ranef = TRUE} a list including "info" and an 
extra component called "ranef" containing the random effects. It is a matrix with the same number of columns
as the data. Each column contains the randome effects of each variable. 
}

\references{
D.C. Montgomery (2001). Design and analysis of experiments (5th Edition). New York: John Wiley & Sons. 

Charles S. Davis (2002). Statistical methods for the analysis of repeated measures. New York: Springer-Verlag.

Demidenko E. (2013). Mixed Models: Thoery and Applications with R 2nd Edition). New Jersey: 
John Wiley & Sons (Excellent book).
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{varcomps.mle}, \link{colrint.regbx}
}
}

\examples{
## example taken from Montgomery, page 514-517.
y <- c(98, 97, 99, 96, 91, 90, 93, 92,
96, 95, 97, 95, 95, 96, 99, 98)
y <- matrix(y)
id <- rep(1:4, each = 4)

x <- rmvnorm(100, numeric(20), diag(rexp(20)) )
id <- rep(1:25, each = 4)
n <- 25  ;  d <- 4
a <- colvarcomps.mom(x, id) 
mean(a[, 4]<0 & a[, 5]>0)  
b <- colvarcomps.mle(x, id) 
x <- NULL
}
