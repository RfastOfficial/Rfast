\name{Polyserial correlation}
\alias{poly.cor}
\title{
Polyserial correlation
}

\description{
Polyserial correlation.
}

\usage{
poly.cor(x, y)
}

\arguments{
\item{x}{
The continuous variable.
}
\item{y}{
The ordinal variable, a numeric vector with numbers starting from 1.
}
}

\details{
The polyserial correlation between a continuous and an ordinal variable is calculated. The function is not 
super fast, yet is faster than other implementations we found. 
}

\value{
A list including:
\item{est}{
A vector with the polyserial correlation and its estimated variance.
}
\item{test}{
A vector with the test statistic and its associated p-value.
}
}

\references{
Olsson U., Drasgow F. and Dorans N. J. (1982). 
The polyserial correlation coefficient. Psychometrika, 47(3):337-347.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{correls}, \link{Table}
}
}

\examples{
x <- rnorm(100)
y <- rpois(100, 10) + 1
res<-poly.cor(x, y)
}


