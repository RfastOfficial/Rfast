\name{Norm of a matrix}
\alias{Norm}
\title{
Norm of a matrix
}

\description{
Norm of a matrix.
}

\usage{
Norm(x, type = "F")
}

\arguments{
\item{x}{
A matrix with numbers.
}
\item{type}{
The type of norm to be calculated. The default is "F" standing for Frobenius norm ("f" in R's norm). The other 
options are "C" standing for the one norm ("o" in R's norm), "R" for the identiy norm ("I" in R's norm) and "M" 
for the maximum modulus among elements of a matrix ("M" in R's norm)
}
}

\value{
A number, the norm of the matrix.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{Dist}, \link{dista}, \link{colmeans}
}
}

\examples{
x <- matrix( rnorm(10 * 10), ncol = 10 )
res<-Norm(x, "F")
res<-norm(x, "F")
res<-Norm(x, "M")
res<-norm(x, "M")
}

\keyword{ Norm of a matrix }
