\name{Column and row-wise ranks}
\alias{colRanks}
\alias{rowRanks}
\title{
Column and row-wise ranks
}

\description{
Column and row-wise ranks.
}

\usage{
colRanks(x,method = "average",descending = FALSE,stable = FALSE, parallel = FALSE)
rowRanks(x,method = "average",descending = FALSE,stable = FALSE, parallel = FALSE)
}

\arguments{
\item{x}{
A mumerical matrix with the data.
}
\item{parallel}{
A boolean value for parallel version.
}
\item{method}{
a character string for choosing method. Must be one of "average", "min", "max", "first".
}
\item{descending}{
A boolean value (TRUE/FALSE) for sorting the vector in descending order. By default sorts the vector in ascending.
}
\item{stable}{
A boolean value (TRUE/FALSE) for choosing a stable sort algorithm. Stable means that discriminates on the same elements. Only for the method "first".
}
}

\details{
For each column or row of a matrix the ranks are calculated and they are returned. The initial matrix is gone. 
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

x<-a1<-a2<-b1<-b2<-NULL
}


