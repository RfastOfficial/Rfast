\name{MLE of the inverted Dirichlet distribution}
\alias{invdir.mle}
\title{
MLE of the inverted Dirichlet distribution
}

\description{
MLE of the inverted Dirichlet distribution.
}

\usage{
invdir.mle(x, tol = 1e-07)
}

\arguments{
\item{x}{
A matrix with strictly positive data (no zeros are allowed).
}
\item{tol}{
The tolerance level up to which the maximisation stops.
}
}

\details{
Maximum likelihood estimation of the parameters of the inverted  is performed via Newton-Raphson. We took the initial values suggested by Bdiri T. and Bouguila N. (2012) and modified them a bit. 
}

\value{
A list including:
\item{iters}{
The number of iterations required by the Newton Raphson.
}
\item{loglik}{
The value of the log-likelihood.
}
\item{param}{
The estimated parameters.
}
}

\references{
Bdiri T. and Bouguila N. (2012). Positive vectors clustering using inverted Dirichlet finite mixture models. Expert Systems with Applications, 39(2): 1869-1882.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.nr2}, \link{multinom.mle}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
for(i in 1:10) invdir.mle(x)
res<-invdir.mle(x)
}

\keyword{ Inverted Dirichlet distribution }
\keyword{ positive multivariate data }
