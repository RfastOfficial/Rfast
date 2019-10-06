\name{Hypothesis test for von Mises-Fisher distribution over Kent distribution}
\alias{fish.kent}
\title{
Hypothesis test for von Mises-Fisher distribution over Kent distribution
}

\description{
The null hypothesis is whether a von Mises-Fisher distribution fits the data well, and the
altenrative is that the Kent distribution is more suitable.
}

\usage{
fish.kent(x, logged = FALSE)
}

\arguments{
\item{x}{
A numeric matrix containing the data as unit vectors in Euclidean coordinates.
}
\item{logged}{
If you want the logarithm of the p-value ot be returned set this to TRUE. 
}
}

\details{
Essentially it is a test of rotational symmetry, whether Kent's ovalness parameter (beta) is equal to zero. 
This works for spherical data only.
}

\value{
A vector with two elements, the value of the test statistic and its associated p-value.
}

\references{
Rivest, L. P. (1986). Modified Kent's statistics for testing goodness of fit for the Fisher distribution in 
small concentrated samples. Statistics & probability letters, 4(1): 1-4.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.mle}, \link{iag.mle}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
res<-fish.kent(x)
x <- NULL
}

\keyword{ Goodness of fit test }
\keyword{ von Mises-Fisher distribution }
\keyword{ Kent distribution }

