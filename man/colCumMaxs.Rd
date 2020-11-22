\name{Colum-wise cumulative operations (sum, prod, min, max)}
\alias{colCumSums}
\alias{colCumProds}
\alias{colCumMins}
\alias{colCumMaxs}

\title{
Colum-wise cumulative operations (sum, prod, min, max)
}

\description{
Colum-wise cumulative operations (sum, prod, min, max).
}

\usage{
colCumSums(x)
colCumProds(x)
colCumMins(x)
colCumMaxs(x)
}

\arguments{
\item{x}{
A numerical matrix.
}
}

\details{
Cumulative mins, maxs, sums and prods are returned. 
}

\value{
A matrix with the results. It has one row less than the initial matrix.
}

%\references{
%
%}

\author{
Manos Papadakis and Michail Tsagris

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> 
and Michail Tsagris <mtsagris@uoc.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{colsums}, \link{colMedians}, \link{colVars}
}
}

\examples{
x <- matrnorm(10, 10)
res<-colCumSums(x)
res<-colCumMins(x)
res<-colCumMaxs(x)
res<-colCumProds(x)
}


