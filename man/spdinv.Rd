\name{Inverse of a symmetric positive definite matrix}
\alias{spdinv}
\title{
Inverse of a symmetric positive definite matrix
}

\description{
Inverse of a symmetric positive definite matrix.
}

\usage{
spdinv(A)
}

\arguments{
\item{A}{
A square positive definite matrix.
}
}

\details{
After calculating the Cholesky decomposition of the matrix we use this upper triangular matrix to 
invert the original matrix.
}

\value{
The inverse of the input matrix.
}

\references{
http://econ.nsysu.edu.tw/ezfiles/124/1124/img/Chapter17_MaximumLikelihoodEstimation.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cholesky}, \link{cova}
}
}

\examples{
s <- cova( as.matrix(iris[, 1:4]) )
res<-spdinv(s)
res<-solve(s)
}

\keyword{ Inverse matrix }
\keyword{ positive definite }
\keyword{ symmetric matrix }