\name{Column and row-wise means of a matrix}
\alias{colmeans}
\alias{colmeans.matrix}
\alias{colmeans.data.frame}
\alias{rowmeans}
\alias{colhameans}
\alias{rowhameans}
\title{
Column and row-wise means of a matrix
}

\description{
Column and row-wise means of a matrix.
}

\usage{
colmeans(x, parallel = FALSE)
\method{colmeans}{matrix}(x, parallel = FALSE)
\method{colmeans}{data.frame}(x, parallel = FALSE)
rowmeans(x)
colhameans(x, parallel = FALSE)
rowhameans(x)
}

\arguments{
\item{x}{
A numerical matrix or data.frame with data.
}
\item{parallel}{
Do you want to do it in parallel in C++? TRUE or FALSE.
}
}


\value{
A vector with the column or row arithmetic or harmonic means.
}

%\references{
%Shabalin, A.A. (2012). Matrix eQTL: Ultra fast eQTL analysis %via large matrix operations. Bioinformatics 28(10): 1353-1358.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colsums}, \link{rowsums}, \link{colMins}, \link{colMedians}, \link{colMads}
}
}

\examples{
x <- matrix(rpois(100 * 100, 10),ncol = 100)
x1 <- colmeans(x)
x2 <- colMeans(x)
all.equal(x1,x2)

x1 <- rowmeans(x)
x2 <- rowMeans(x)
all.equal(x1,x2)
system.time( colhameans(x) )
system.time( rowhameans(x) ) 

x<-x1<-x2<-NULL
}

\keyword{ Column means }
\keyword{ row means }
\keyword{ harmonic means }