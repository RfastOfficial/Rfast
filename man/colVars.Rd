\name{Column and row-wise variances and standard deviations}
\alias{colVars}
\alias{rowVars}
\title{
Column and row-wise variances and standard deviations of a matrix
}

\description{
Column and row-wise variances and standard deviations of a matrix
}

\usage{
colVars(x, std = FALSE, na.rm = FALSE, parallel = FALSE)
rowVars(x, std = FALSE, na.rm = FALSE, parallel = FALSE)
}

\arguments{
\item{x}{
A matrix with the data.
}
\item{std}{
A boolean variable specyfying whether you want the variances (FALSE) 
or the standard deviations (TRUE) of each column.
}
\item{na.rm}{
TRUE or FAlSE for remove NAs if exists.
}
\item{parallel}{
Should parallel implentations take place in C++? The default value is FALSE.
}
}

\details{
We found this on stackoverflow which was created by David Arenburg. We then modified the function to match the sums type formula of the variance, which is faster.
}

\value{
A vector with the column variances or standard deviations.
}

\author{
Michail Tsagris and Manos Papadakis.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{colmeans}, \link{colMedians}, \link{colrange}
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100 )
a2 <- colVars(x)
x<-a2<-NULL
}

\keyword{ Column-wise variances }
\keyword{ row-wise variances }
\keyword{ variances of many samples }
