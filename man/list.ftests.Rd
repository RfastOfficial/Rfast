\name{Many F-tests with really huge matrices}
\alias{list.ftests}
\title{
Many F-tests with really huge matrices
}

\description{
Many F-tests with really huge matrices.
}

\usage{
list.ftests(x, logged = FALSE)
}

\arguments{
\item{x}{
A list with many big size matrices. Each element of the list contains a matrix. This is the \code{\link{ftests}}
function but with really huge matrices, which cannot be loaded into R as a single matrix.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
The Welch's F-test (without assuming equal variances) is performed just like in the "ftests" function.
The difference is that you have a really huge matrix which you cannot load into R. In the "ftests" function,
the argument "ina" denotes the different groups. Here, you "cut" the matrix into smaller ones, each of which
denotes a different group and put them in a list.
}

\value{
A matrix with the test statistic and the p-value of each test.
}

\references{
B.L. Welch (1951). On the comparison of several mean values: an alternative approach. Biometrika, 38(3/4), 330-336.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ftests}, \link{ttests}
}
}

\examples{
x <- matrnorm(300, 500)
ina <- rbinom(300, 2, 0.6) + 1
a <- list()
a[[ 1 ]] <- x[ina == 1, ]
a[[ 2 ]] <- x[ina == 2, ]
a[[ 3 ]] <- x[ina == 3, ]
mod <- list.ftests(a)
z <- NULL
a <- NULL
}

\keyword{ F-tests }
\keyword{ analysis of variance }
\keyword{ huge datasets }
