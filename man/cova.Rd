\name{Covariance and correlation matrix}
\alias{cova}
\alias{cora}
\title{
Fast covariance and correlation matrix calculation
}

\description{
Fast covariance and correlation matrix calculation.
}

\usage{
cova(x, center = FALSE, large = FALSE)

cora(x, large = FALSE)
}

\arguments{
\item{x}{
A matrix with data. It has to be matrix, if it is data.frame for example the function does not turn it into a matrix.
}
\item{center}{
If you want to center the data prior to applying the cross product of the mateix set this equal to TRUE, otherwise leave it NULL.
}
\item{large}{
If you have large matrices, with thousands of rows and or many tens or hundreds of columns set this equal to TRUE in 
order to use Rfast's \code{\link{Crossprod}} or \code{\link{\Tcrossprod}} functions. These functions are twice or up to
3 times faster than the correpsonding built-in functions. 
}
}

\details{
The calculations take place faster than the built-in functions \code{\link{cor}} 
as the number of variables increases. This is true if the number of variables is high, 
say from 500 and above. The "cova" on the other hand is always faster. For the "cova" in specific, 
we have an option to center the data prior to the cross product. This can be more stable if you 
have many tens of thousands of rows due to numerical issues that can arise. 

For the correlation matrix we took the code from here

https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
}

\value{
The covariance or the correlation matrix.
}

\author{
Michail Tsagris and Manos Papadakis <papadakm95@gmail.com>.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{colVars}, \link{cor}, \link{cov}
}
}

\examples{
x <- matrnorm(100, 40)
s1 <- cov(x) 
s2 <- cova(x)
all.equal(s1, s2)
x <- NULL
}

\keyword{ Covariance matrix }
