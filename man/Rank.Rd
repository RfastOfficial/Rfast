\name{Column and row-wise ranks}
\alias{colRanks}
\alias{rowRanks}
\alias{Rank}
\title{
Column and row-wise ranks
}

\description{
Column and row-wise ranks.
}

\usage{
colRanks(x,method = "average",descending = FALSE,stable = FALSE, parallel = FALSE, cores = 0)
rowRanks(x,method = "average",descending = FALSE,stable = FALSE, parallel = FALSE, cores = 0)
Rank(x,method = "average",descending = FALSE,stable = FALSE, parallel = FALSE)
}

\arguments{
\item{x}{
A mumerical matrix or data.frame with the data.
}
\item{parallel}{
A boolean value for parallel version.
}
\item{method}{
a character string for choosing method. Must be one of:
\itemize{
\item average : a permutation with theri mean values at each index set of ties
\item min : a permutation with minimum values at each index set of ties
\item max : a permutation with maximum values at each index set of ties
\item first : a permutation with increasing values at each index set of ties
\item random : a permutation with random values at each index set of ties
}
}
\item{descending}{
A boolean value (TRUE/FALSE) for sorting the vector in descending order. By default sorts the vector in ascending.
}
\item{stable}{
A boolean value (TRUE/FALSE) for choosing a stable sort algorithm. Stable means that discriminates on the same elements. Only for the method "first".
}
\item{cores}{
Number of cores to use for parallelism. Valid only when argument parallel is set to TRUE. 
Default value is 0 and it means the maximum supported cores.
}
}

\details{
For each column or row of a matrix the ranks are calculated and they are returned. 
}

\value{
A matrix with the column or row-wise ranks.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{Rank}, \link{correls}
}
}

\examples{
x <- matrnorm(100, 10)
a1 <- colRanks(x)
a2 <- apply(x, 2, rank)
b1 <- rowRanks(x)
b2 <- apply(x, 1, rank)

a1 <- Rank(x[,1])
a1 <- rank(x[,1])

x<-a1<-a2<-b1<-b2<-NULL
}


