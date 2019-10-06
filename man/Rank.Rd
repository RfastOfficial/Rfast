\name{Ranks of the values of a vector}
\alias{Rank}
\title{
Ranks of the values of a vector
}

\description{
Ranks of the values of a vector.
}

\usage{
Rank(x,method = "average",descending = FALSE)
}

\arguments{
\item{x}{
A numerical vector with data.
}
\item{method}{
a character string for choosing method. Must be one of "average", "min", "max".
}
\item{descending}{
A boolean value (TRUE/FALSE) for sorting the vector in descending order. By default sorts the vector in ascending.
}
}

\details{
The ranks of the values are returned, the same job as "rank". If you want you can choose descending/ascending order for all methods.
}

\value{
A vector with the ranks of the values.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{colRanks}, \link{correls}
}
}

\examples{
x <- rnorm(100)
a1 <- Rank(x)
a2 <- rank(x)
}


