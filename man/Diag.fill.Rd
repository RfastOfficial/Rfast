\name{Diagonal Matrix}
\alias{Diag.fill}
\alias{Diag.matrix}
\title{
Diagonal Matrix
}

\description{
Fill the diagonal of a matrix or create a diagonal and initialize it with a specific value.
}

\usage{
Diag.fill(x,v=0)
Diag.matrix(len,v=0)
}

\arguments{
\item{x}{
A matrix with data.
}
\item{len}{
Number of columns or rows.
}
\item{v}{
Value or vector to initialize the diagonal of a matrix.By default "v=0".
}
}


\value{
Diag.fill returns a diagonal matrix where all the elements in the diagonal are equal to "v". 

Diag.matrix returns a diagonal matrix where has dimension "len,len" and all the elements in the diagonal are equal to "v". It is fast for huge matrices with dimensions more than [row,col] = [500,500]
}

%\references{
%Shabalin, A.A. (2012). Matrix eQTL: Ultra fast eQTL analysis %via large matrix operations. Bioinformatics 28(10): 1353-1358.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{rowMins}, \link{colFalse}, \link{nth}, \link{rowrange}, \link{rowMedians}, \link{rowVars}, \link{colSort}, \link{rowSort}, \link{colTrue}
}
}

\examples{
x <- matrix(rbinom(100*100,1,0.5),100,100)

f <- Diag.fill(x,1)
f <- Diag.fill(x,1:100) ##equals to diag(x)<-1:100
f <- Diag.matrix(100,1) ##equals to diag(1,100,100)
f <- Diag.matrix(100,1:100) ##equals to diag(1:100,100,100)

f<-x<-NULL
}

\keyword{ Diagonal Matrix }
\keyword{ Create - Fill }