\name{Natural logarithm of the gamma function and its derivatives}
\alias{Lgamma}
\alias{Digamma}
\alias{Trigamma}
\title{
Natural logarithm of the gamma function and its derivatives.
}

\description{
Natural logarithm of the gamma function and its derivatives.
}

\usage{
Lgamma(x)
Digamma(x)
Trigamma(x)
}

\arguments{
\item{x}{
A numerical matrix or vector with positive numbers in either case.
}
}

\details{
We have spotted that the time savings come when there are more than 50 elements, with vector or matrix.
}

\value{
The matrix or the vector with the resulting values.
}

\references{
Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical Functions. New York: Dover. https://en.wikipedia.org/wiki/Abramowitz_and_Stegun provides links to the full text which is in public domain.
Chapter 6: Gamma and Related Functions. 
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{beta.mle}, \link{diri.nr2}
}
}

\examples{
x <- matrix( rnorm(500 * 500), ncol = 500 )
 a1 <- Lgamma(x) 
a2 <- lgamma(x) 
all.equal(as.vector(a1), as.vector(a2))

a1 <- Digamma(x) 
a2 <- digamma(x) 
all.equal(as.vector(a1), as.vector(a2))

x<-a1<-a2<-NULL
}

\keyword{ Logarithm of gamma function }
\keyword{ derivatives }
\keyword{ digamma function }
\keyword{ trigamma function }

