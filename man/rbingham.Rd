\name{Simulation of random values from a Bingham distribution with any symmetric matrix}
\alias{rbingham}
\title{
Simulation of random values from a Bingham distribution with any symmetric matrix
}
\description{
Simulation of random values from a Bingham distribution with any symmetric matrix.
}

\usage{
rbingham(n, A)
}

\arguments{
\item{n}{
Sample size.
}
\item{A}{
A symmetric matrix. 
}
}

\details{
The eigenvalues of the q x q symmetric matrix A are calculated and the smallest of 
them is subtracted from the rest. The q - 1 non zero eiqenvalues are then passed to 
\code{\link{rbing}}. The generated data are then right multiplied by \eqn{V^T}, where
V is the matrix of eigenvectors of the matrix A.
}

\value{
A matrix with the simulated data.
}

\references{
Kent J.T., Ganeiber A.M. and Mardia K.V. (2013).
A new method to simulate the Bingham and related distributions
in directional data analysis with applications.
http://arxiv.org/pdf/1310.8110v1.pdf

C.J. Fallaize and T. Kypraios (2014).
Exact Bayesian Inference for the Bingham Distribution.
Statistics and Computing (No volum assigned yet).
http://arxiv.org/pdf/1401.2894v1.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and 
Manos Papadakis <papadakm95@gmail.com>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{rvmf}}
}

\examples{
A <- cov( iris[, 1:4] )
x <- rbingham(100, A)
x
}

