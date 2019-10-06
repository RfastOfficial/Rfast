\name{Random values simulation from a von Mises distribution}
\alias{rvonmises}
\title{
Random values simulation from a von Mises distribution
}

\description{
It generates random vectors following the von Mises distribution. The data can be spherical or hyper-spherical.
}

\usage{
rvonmises(n, m, k, rads = TRUE)
}

\arguments{
\item{n}{
The sample size.
}
\item{m}{
The mean angle expressed in radians or degrees.
}
\item{k}{
The concentration parameter. If k is zero the sample will be generated from the uniform distribution over \eqn{(0, 2\pi)}.
}
\item{rads}{
If the mean angle is expressed in radians, this should be TRUE and FALSE otherwise. The simulated data will be expressed in radians or degrees depending on what the mean angle is expressed.
}
}

\details{
The mean direction is transformed to the Euclidean coordinates (i.e. unit vector) and then the fvmf function is employed. It uses a rejection smapling as suggested by Andrew Wood in 1994. I have mentioned the description of the algorithm as I found it in Dhillon and Sra in 2003. Finally, the data are transformed to radians or degrees.
}

\value{
A vector with the simulated data.
}

\references{
Wood, A. T. (1994). Simulation of the von Mises Fisher distribution. Communications in statistics-simulation 
and computation, 23(1): 157-164. 

Dhillon, I. S., & Sra, S. (2003). Modeling data using directional distributions. Technical Report TR-03-06, 
Department of Computer Sciences, The University of Texas at Austin.
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.75.4122&rep=rep1&type=pdf
}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm85@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vm.mle}, \link{rvmf}
}
}

\examples{
x <- rvonmises(1000, 2, 25, rads = TRUE)
res<-vm.mle(x)
}

\keyword{ Random values simulation }
\keyword{ von Mises-Fisher distribution }

