\name{All k possible combinations from n elements}
\alias{comb_n}
\title{
All k possible combinations from n elements
}

\description{
All k possible combinations from n elements.
}

\usage{
comb_n(n, k,simplify=TRUE)
}

\arguments{
\item{n}{
A positive \bold{INTEGER} number or a vector with numbers. 
}
\item{k}{
A positive integer number at most equal to n or at most equal to the length of n, if n is a vector.
}
\item{simplify}{
A logical value for return List instead of matrix.
}
}


\value{
A matrix with k columns and rows equal to the number of possible unique combinations of n with k elements.
If simplify is set to TRUE then a list with k values where each value has length equal to the number of 
possible unique combinations of n with k elements.
}

\references{
Nijenhuis A. and Wilf H.S. (1978). Combinatorial Algorithms for Computers and Calculators. Academic Press, NY. 
}

\author{
Manos Papadakis and Marios Dimitriadis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> and Marios Dimitriadis <kmdimitriadis@gmail.com>.
}


\seealso{
\code{\link{nth}, \link{colMaxs}, \link{colMins}, \link{colrange}
}
}

\examples{
comb_n(20, 4)
combn(20, 4)
x <- rnorm(5)
res<-comb_n(x, 3)
}

\keyword{ Combinatorics }
\keyword{ All possibe combinations }