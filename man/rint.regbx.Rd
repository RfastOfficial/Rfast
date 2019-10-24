\name{Random intercepts linear mixed models}
\alias{rint.reg}
\alias{rint.regbx}
\title{
Random intercepts linear mixed models
}

\description{
Random intercepts linear mixed models (for balanced data with a single identical covariate).
}

\usage{
rint.reg(y, x, id ,tol = 1e-08, ranef = FALSE, maxiters = 100)
rint.regbx(y, x, id)  
}

\arguments{
\item{y}{
A numerical vector with the data. The subject values.
}
\item{x}{
For the case of "rint.reg" this can be a vector or a numerical matrix with data. In the case of "rint.regbx" this is a numerical 
vector with the same length as y indicating the fixed predictor variable. Its values are the same for all levels of y. 
An example of this x is time which is the same for all subjects.
}
\item{id}{
A numerical variable with 1, 2, ... indicating the subject. 
}
\item{tol}{
The tolerance level to terminate the generalised elast squares algorithm. 
}
\item{ranef}{
If you want to obtain the random effects (random intercepts) set this equal to TRUE.
}
\item{maxiters}{
The max number of iterations that can take place in a regression.
}
}

\details{
Random intercepts linear mixed models with compound covariance structure is fitted in both functions. The "rint.reg" allows any numerical 
matrix, with balanced or unbalanced data. See Demidenko (2013, pg. 65-67) for more information. 

The "rint.regbx" is a special case of a balanced random intercepts model with a compound symmetric covariance matrix and one single covariate 
which is constant for all replicates. An example, is time, which is the same for all subjects. Maximum likelihood estimation has been performed. 
In this case the mathematics exist in a closed formula (Demidenko, 2013, pg. 67-69). 
}

\value{
A list including:
\item{info}{
A vector with the random intercepts variance (between), the variance of the errors (within), the log-likelihood, the deviance 
(twice the log-likelihood) and the BIC. In the case of "rint.reg" it also includes the number of iterations required by the 
generalised least squares. 
}
\item{be}{
The estimated regression coefficients, which in the case of "rint.regbx" are simply two: the constant and the slope (time effect). 
}
\item{ranef}{
The random intercepts effects. 
}
}

\references{
Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons (excellent book).
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{rm.lines}, \link{varcomps.mom}, \link{colvarcomps.mom}
}
}

\examples{
\dontrun{
y <- rnorm(100) 
x <- rnorm(10)
x <- rep(x, 10)
id <- rep(1:10, each = 10)
system.time( for (i in 1:40) a <- rint.reg(y, x, id) )
}
}

\keyword{ Linear mixed models }
\keyword{ balanced design }
\keyword{ maximum likelihood estimation }
