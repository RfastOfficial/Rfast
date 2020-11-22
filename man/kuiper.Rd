\name{Uniformity test for circular data}
\alias{kuiper}
\alias{watson}
\title{
Uniformity tests for circular data
}

\description{
Hypothesis tests of uniformity for circular data.
}

\usage{
kuiper(u)

watson(u)
}

\arguments{
\item{u}{
A numeric vector containing the circular data which are expressed in radians.
}
}

\details{
These tests are used to test the hypothesis that the data come from a circular uniform distribution. 
}

\value{
A vector with two elements, the value of the test statistic and its associated p-value.
}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, pg. 153-55 (Kuiper's test) & 
156-157 (Watson's test).
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.mle}, \link{rvonmises}
}
}

\examples{
x <- rvonmises(n = 50, m = 2, k = 10)
res<-kuiper(x)
res<-watson(x)
x <- runif(50, 0, 2 * pi)
res<-kuiper(x)
res<-watson(x)
}

\keyword{ Circular data}
\keyword{ uniformity test }
\keyword{ hypothesis testing }



