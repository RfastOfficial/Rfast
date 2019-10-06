\name{ Deep copy and printing of an environment }
\alias{ env.copy }
\alias{ print.environment }
\title{
Deep copy and printing of an environment
}

\description{
Deep copy and printing of an environment.
}

\usage{
env.copy(x,all.names=FALSE)
\method{print}{environment}(x,all.names=FALSE,...)
}

\arguments{
\item{x}{
An environment object.
}
\item{all.names}{
An logical value (TRUE or FALSE):

env.copy: copy all the hidden variables or not
print.environment: print all the hidden variables or not
}
\item{...}{
Anything the user want.
}
}

\details{
env.copy: deep copy of the environment object.
print.environment: print the name, type and length or dimension of each variable inside environment. For printing the hidden objects set second argument to TRUE, by default it is FALSE.
}

\value{
env.copy: A copy of the first argument.
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colShuffle}, \link{colVars}, \link{colmeans}, \link{read.directory}
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
\keyword{ Printing }
