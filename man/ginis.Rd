\name{Many Gini coefficients}
\alias{ginis}
\title{
Many Gini coefficients
}

\description{
Many Gini coefficients.
}

\usage{
ginis(x)
}

\arguments{
\item{x}{
A matrix with non negative data. The rows are observations and the columns denote the variables.
}
}

\details{
We have implemented the fast version of the Gini coefficient. See \href{https://en.wikipedia.org/wiki/Gini_coefficient}{wikipedia} for more details.
}

\value{
A vector with the Gini coefficient, one for each variable.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colskewness}, \link{colmeans}, \link{corpairs}
}
}

\examples{
x <- matrix( rpois(500 * 1000, 1000), ncol = 1000 )
a <- ginis(x)
}

\keyword{ Gini coefficient }


