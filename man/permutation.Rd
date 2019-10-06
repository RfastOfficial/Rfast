\name{ Permutation }
\alias{permutation}
\alias{permutation.next}
\alias{permutation.prev}
\alias{bincomb}
\title{
Permutation
}

\description{
Permute the given vector.
}

\usage{
permutation(x, nperm = gamma(length(x)+1))
permutation.next(x, nperm = gamma(length(x)+1))
permutation.prev(x, nperm = gamma(length(x)+1))
bincomb(n)
}

\arguments{
\item{x}{
A numeric vector with data. 
}
\item{nperm}{
An integer value for returning specific number of combinations. By defualt is set to all combinations. Must be between \bold{0<=nperm<=gamma(length(x)+1)}
}
\item{n}{
An integer value for the length of the binary number.
}
}

\details{
This function implements "Permutation", which means all the possible combinations. In the permutation.next and permutation.prev
if there aren't possible combinations it returns the same vector. "Binary Combinations" for "bincomb", means all the possible combinations for the binary number with length "n".
}

\value{
Returns a matrix with all possible combinations of the given vector or a matrix row with one possible combinations.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>
}

\seealso{
\code{\link{combn},\link{comb_n}
}
}

\examples{
y <- rnorm(3)
b <- permutation(y)
b <- permutation.next(y)
b <- permutation.prev(y)
g <- bincomb(3)
}

\keyword{ Permutation Function }