\name{Column and row-wise range of values of a matrix}
\alias{colrange}
\alias{rowrange}
\title{
Column and row-wise range of values of a matrix
}

\description{
Column and row-wise range of values of a matrix.
}

\usage{
colrange(x, cont = TRUE, parallel = FALSE,cores = 0)
rowrange(x, cont = TRUE)
}

\arguments{
\item{x}{
A numerical matrix or data.frame with data.
}
\item{parallel}{
Execute algorithm in parallel for data.frame.
}
\item{cont}{
If the data are continuous, leave this TRUE and it will return the range of values for each variable (column). If the data are integers, categorical, or if you want 
to find out the number of unique numbers in each column set this to FALSE.
}
\item{cores}{
Number of cores to use for parallelism. Valid only when argument parallel is set to TRUE. 
Default value is 0 and it means the maximum supported cores.
}
}


\value{
A vector with the relevant values.
}

%\references{
%Shabalin, A.A. (2012). Matrix eQTL: Ultra fast eQTL analysis %via large matrix operations. Bioinformatics 28(10): 1353-1358.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colMins}, \link{colMaxs}, \link{rowMins}, \link{rowMaxs}, \link{nth}, \link{colMedians}, \link{colVars}, \link{colSort}, \link{rowSort}
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100 )

a1 <- colrange(x) 
a2 <- apply(x, 2, function(x) diff( range(x)) ) 
all.equal(a1, a2)

a1 <- rowrange(x) 
a2 <- apply(x, 1, function(x) diff( range(x)) ) 
all.equal(a1, a2)

x<-a1<-a2<-NULL
}

\keyword{ Column-wise ranges}

