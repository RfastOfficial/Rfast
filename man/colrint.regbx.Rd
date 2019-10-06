\name{Many random intercepts LMMs for balanced data with a single identical covariate.}
\alias{colrint.regbx}
\title{
Many random intercepts LMMs for balanced data with a single identical covariate
}

\description{
Many random intercepts LMMs for balanced data with a single identical covariate.
}

\usage{
colrint.regbx(y, x, id)  
}

\arguments{
\item{y}{
A numerical matrix with the data. The subject values.
}
\item{x}{
A numerical vector with the same length as the number of rows of y indicating the fixed predictor variable. 
Its values are the same for all levels of y. An example of this x is time which is the same for all subjects.
}
\item{id}{
A numerical variable with 1, 2, ... indicating the subject. 
}
}

\details{
This is a special case of a balanced random intercepts model with a compound symmetric covariance matrix and one 
single covariate which is constant for all replicates. An example, is time, which is the same for all subjects. 
Maximum likelihood estimation has been performed. In this case the mathematics exist in a closed formula 
(Demidenko, 2013, pg. 67-69). 

This is the generalistion of \code{\link{rint.regbx}} to matrices. Assume you have many observations, gene 
expressions over time for example, and you want to calculate the random effects or something else for each 
expression. Instead of using a "for" loop with \code{\link{rint.regbx}} function we have used amtrix operations to 
make it even faster.
}

\value{
A list including:
\item{info}{
A matrix with the random intercepts variance (between), the variance of the errors (within), the log-likelihood, the 
deviance (twice the log-likelihood) and the BIC. In the case of "rint.reg" it also includes the number of iterations 
required by the generalised least squares. 
}
\item{be}{
The estimated regression coefficients, which in the case of "rint.regbx" are simply two: the constant and the slope 
(time effect). 
}
\item{ranef}{
A matrix with random intercepts effects. Each row corresponds to a column in y. Instead of having a matrix with the 
same number of columns as y we return a transposed matrix. 
}
}

\references{
Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons 
(excellent book).
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{colvarcomps.mle}, \link{rint.regbx}, \link{rm.lines}, \link{varcomps.mom}, \link{rint.reg}
}
}

\examples{
y <- matrix( rnorm(100 * 50), ncol = 50)
id <- rep(1:20, each = 5)
x <- rep(1:10, 10)
system.time( a<- colrint.regbx(y, x, id) )
x <- NULL
}

\keyword{ Linear mixed models }
\keyword{ balanced design }
\keyword{ maximum likelihood estimation }
