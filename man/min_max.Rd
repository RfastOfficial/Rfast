\name{minimum and maximum }
\alias{min_max}
\title{
Minimum and maximum of a vector
}

\description{
Minimum and maximum of a vector.
}

\usage{
min_max(x,index=FALSE, percent = FALSE)
}

\arguments{
\item{x}{
A numerical vector with data. NAs are handled naturally.
}
\item{index}{
A boolean value for the indices of the minimum and the maximum value.
}
\item{percent}{
A boolean value for the percent of the positive and negative numbers.
}
}

%\details{
%
%}

\value{
A vector with the relevant values, min and max.
}

%\references{
%Shabalin, A.A. (2012). Matrix eQTL: Ultra fast eQTL analysis %via large matrix operations. Bioinformatics 28(10): 1353-1358.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{rowMins}, \link{rowMaxs}, \link{nth}, \link{colrange}, \link{colMedians},\link{colSort}, \link{rowSort}
}
}

\examples{
x <- rnorm(100 * 500)
s1 <- min_max(x) 
s2 <- c(min(x), max(x)) 

}

\keyword{ minimum }
\keyword{ maximum }

