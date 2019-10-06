\name{Natural logarithm of the beta function}
\alias{Lbeta}
\title{
Natural logarithm of the beta function
}

\description{
Natural logarithm of the beta function.
}

\usage{
Lbeta(x, y)
}

\arguments{
\item{x}{
A numerical matrix, or a vector or just a number with positive numbers in either case.
}
\item{y}{
A numerical matrix, or a vector or just a number with positive numbers in either case. The dimensions of y 
must match those of x.
}
}

\details{
The function is faster than R's lbeta when the dimensions of x any are large. If you have only two numbers, 
then lbeta is faster. But if you have for example two vectors of 1000 values each, Lbeta becomes two times 
faster than lbeta.
}

\value{
The matrix, vector or number with the resulting values.
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
\code{\link{Lgamma}, \link{beta.mle}, \link{diri.nr2}
}
}

\examples{
x <- rexp(1000)
y <- rexp(1000)
a1 <- Lbeta(x, y)

x<-y<-a1<-NULL
}

\keyword{ Beta function }
\keyword{ logarithm }

