\name{Simulation of random values from a Bingham distribution}
\alias{rbing}
\title{
Simulating from a Bingham distribution}
\description{
Simulation from a Bingham distribution using the code suggested by Kent et al. (2013).
}

\usage{
rbing(n, lam)
}

\arguments{
\item{n}{
Sample size.
}
\item{lam}{
Eigenvalues of the diagonal symmetric matrix of the Bingham distribution. 
See details for more information on this.
}
}

\details{
The user must have calculated the eigenvalues of the diagonal symmetric matrix of the 
Bingham distribution. The function accepts the q-1 eigenvalues only. This means, that 
the user must have subtracted the lowest eigenvalue from the rest and give the non 
zero ones. The function uses rejection sampling. 
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

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and 
Manos Papadakis <papadakm95@gmail.com>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{rvmf}}
}

\examples{
x <- rbing( 100, c(1, 0.6, 0.1) )
x
}

