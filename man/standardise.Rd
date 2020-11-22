\name{Standardisation}
\alias{standardise}
\title{
Standardisation
}

\description{
Standardisation.
}

\usage{
standardise(x, center = TRUE, scale = TRUE)
}

\arguments{
\item{x}{
A matrix with data. It has to be matrix, if it is data.frame for example the function does not turn it into a matrix.
}
\item{center}{
Should the data be centred as well? TRUE or FALSE.
}
\item{scale}{
Should the columns have unit variance, yes (TRUE) or no (FALSE)?
}
}

\details{
Similar to R's built in functions "scale" there is the option for centering or scaling only or both (default).
}

\value{
A matrix with the standardised data.
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{colVars}, \link{colmeans}, \link{colMads}
}
}

\examples{
x <- matrnorm( 100, 100 )
a1 <- scale(x)[1:100, ] 
a2 <- standardise(x) 
all.equal(as.vector(a1), as.vector(a2))
x <- NULL
}

\keyword{ Standardisation }
