\name{ Deep copy }
\alias{ env.copy }
\title{
Deep copy
}

\description{
Deep copy.
}

\usage{
env.copy(x,all.names=FALSE)
}

\arguments{
\item{x}{
An environment object.
}
\item{all.names}{
An logical value (TRUE or FALSE). Copy all the hidden variables or not.
}
}

\details{
Deep copy of the environment object.
}

\value{
A copy of the first argument.
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colShuffle}, \link{colVars}, \link{colmeans}
}
}

\examples{
x <- new.env()
x$imaginary <- NULL
x$real <- NULL

# you can library the package and just press x and R will understand 
# and search automatically for a function to print the environment
x

y <- env.copy(x)

x$real <- 10

x$real == y$real # FALSE

}

\keyword{ Environment }
\keyword{ Deep copy }
