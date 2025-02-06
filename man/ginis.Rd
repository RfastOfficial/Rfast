\name{One and many Gini coefficients}
\alias{gini}
\alias{ginis}
\title{
One and many Gini coefficients
}

\description{
Many Gini coefficients.
}

\usage{
gini(x)
ginis(x)
}

\arguments{
\item{x}{
Either a vector or a matrix with non negative data. The rows are observations and the columns denote the variables.
}
}

\details{
We have implemented the fast version of the Gini coefficient. See \href{https://en.wikipedia.org/wiki/Gini_coefficient}{wikipedia} for more details.
}

\value{
A value with the Gini coefficient or a vector with the Gini coefficient, one for each variable.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colskewness}, \link{colmeans}, \link{corpairs}
}
}

\examples{
x <- matrix( rpois(300 * 100, 100), ncol = 100 )
a <- ginis(x)
}

\keyword{ Gini coefficient }


