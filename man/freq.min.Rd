\name{Minimum and maximum frequencies }
\alias{freq.min}
\alias{freq.max}
\title{
Minimum and maximum frequencies of a vector
}

\description{
Minimum and maximum frequencies of a vector.
}

\usage{
freq.min(x,na.rm = FALSE)
freq.max(x,na.rm = FALSE)
}

\arguments{
\item{x}{
A numerical/integer vector with data but without NAs.
}
\item{na.rm}{
TRUE or FAlSE for remove NAs if exists.
}
}


\value{
A vector with 2 values, the value with minimum/maximum frequency and the frequency.
}

\details{
Those functions are the same with max(table(x) or min(table(x)) but with one exception.
freq.min and freq.max will return also which value has the minimum/maximum frequency.
More Efficient than max(table(x) or min(table(x)).
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> and Marios Dimitriadis <kmdimitriadis@gmail.com>.
}

\seealso{
\code{\link{rowMins}, \link{rowMaxs}, \link{nth}, \link{colrange}, \link{colMedians},\link{colSort}, \link{rowSort}
}
}

\examples{
x <- rnorm(100)
f1 <- freq.min(x) 
f2 <- freq.max(x) 
# f1r <- min(table(x))
# f2r <- max(table(x))
# f1[2]==f1r ## the frequencies are the same
# f2[2]==f2r ## the frequencies are the same

}

\keyword{ minimum frequency}
\keyword{ maximum frequency}