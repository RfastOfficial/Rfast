\name{Partial distance correlation}
\alias{pdcor}
\title{
Partial distance correlation
}

\description{
Partial distance correlation.
}

\usage{
pdcor(x, y, z)
}

\arguments{
\item{x}{
A numerical matrix.
}
\item{y}{
A numerical matrix.
}
\item{z}{
A numerical matrix.
}
}

\details{
The partial distance correlation between x and y conditioning on z is computed.
}

\value{
The partial distance correlation
}

\references{
Szekely G.J. and Rizzo M.L. (2014). Partial Distance Correlation with Methods for Dissimilarities.
Annals of Statistics, 42(6): 2382--2412.

Szekely G. J. and Rizzo M. L. (2023). The Energy of Data and Distance Correlation. 
Chapman and Hall/CRC.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}

\seealso{
\code{ \link{dcov}, \link{dcor.ttest}
}
}

\examples{
x <- as.matrix(iris[1:50, 1:4])
y <- as.matrix(iris[51:100, 1:4])
z <- as.matrix(iris[101:150, 1:4])
#pdcor(x, y, z)
x <- y <- NULL
}

