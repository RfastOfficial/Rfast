\name{Mean - Median absolute deviation of a vector}
\alias{mad2}
\alias{Mad}
\title{
Mean - Median absolute deviation of a vector
}

\description{
Mean - Median absolute deviation of a vector.
}

\usage{
mad2(x,method = "median",na.rm=FALSE)
Mad(x,method = "median",na.rm=FALSE)
}

\arguments{
\item{x}{
A numerical vector. 
}
\item{method}{
A character vector with values "median", for median absolute deviation or "mean", for mean absolute deviation. 
}
\item{na.rm}{
A logical value TRUE/FALSE to remove NAs.
}
}

\value{
The mean absolute deviation. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{colMads}, \link{Median}, \link{colMedians}
}
}

\examples{
x <- Rnorm(1000)
Mad(x)
mad(x)
}




