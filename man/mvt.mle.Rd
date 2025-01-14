\name{MLE of the multivariate t distribution}
\alias{mvt.mle}
\title{
MLE of the multivariate t distribution
}

\description{
MLE of the multivariate t distribution.
}

\usage{
mvt.mle(x, v = 5, tol = 1e-07)
}

\arguments{
\item{x}{
A matrix with numerical data.
}
\item{v}{
The degrees of freedom. Must be a positive number, greater than zero.
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
}

\details{
The location vector, scatter matrix and the value of the log-likelihood is calculated.
}

\value{
A list including:
\item{iters}{The number of iterations required for the EM algorihm to converge.
}
\item{loglik}{The value of the maximised log-likelihood.
}
\item{location}{
The location vector.
}
\item{scatter}{
The scatter matrix.
}
}

\references{
Nadarajah S. and Kotz S. (2008). Estimation methods for the multivariate t distribution.
Acta Applicandae Mathematicae, 102(1):99-118.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{mvnorm.mle}, \link{dmvnorm}, \link{gaussian.nb}
}
}

\examples{
x <- matrnorm(100, 4)
#res<-mvnorm.mle(x)
#res<-mvt.mle(x, v = 5)
#res<-mvt.mle(x, v = 100)
}

