\name{Circular-linear correlation}
\alias{circlin.cor}
\title{
Circular-linear correlation
}
\description{
It calculates the squared correlation between a circular and one or more linear variables.
}
\usage{
circlin.cor(theta, x)
}
\arguments{
\item{theta}{
A circular variable expressed in radians.
}
\item{x}{
The linear variable or a matrix containing many linear variables.
}
}

\details{
The squared correlation between a circular and one or more linear variables is calculated.
}
\value{
A matrix with as many rows as linear variables including:
\item{R-squared}{
The value of the squared correlation.
}
\item{p-value}{
The p-value of the zero correlation hypothesis testing.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{spml.reg}
}
}

\examples{
#phi <- rvonmises(50, 2, 20, rads = TRUE)
#x <- 2 * phi + rnorm(50)
#y <- matrix(rnorm(50 * 5), ncol = 5)
#res<-circlin.cor(phi, x)
#res<-circlin.cor(phi, y)
y <- NULL
}
\keyword{ Circular-linear correlation }
