\name{Column-wise uniformity Watson test for circular data}
\alias{colwatsons}
\title{
Column-wise uniformity tests for circular data
}

\description{
Column-wise uniformity tests for circular data.
}

\usage{
colwatsons(u)
}

\arguments{
\item{u}{
A numeric matrix containing the circular data which are expressed in radians. Each column is a different sample.
}
}

\details{
These tests are used to test the hypothesis that the data come from a circular uniform distribution. 
The Kuiper test is much more time consuming and this is why it not implemented yet. Once we figure out a way to
make it fast, we will incldue it. 
}

\value{
A matrix with two columns, the value of the test statistic and its associated p-value.
}


\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, pg. 153-55 (Kuiper's test) & 
156-157 (Watson's test).
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{watson}, \link{vmf.mle}, \link{rvonmises}
}
}

\examples{
x <- matrix( rvonmises(n = 50 * 10, m = 2, k = 0), ncol = 10 )
res<-colwatsons(x)
x <- NULL
}

\keyword{ Circular data}
\keyword{ uniformity tests }
\keyword{ hypothesis testing }



