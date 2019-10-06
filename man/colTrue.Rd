\name{Column-wise true/false value }
\alias{colTrue}
\alias{colFalse}
\alias{colTrueFalse}
\title{
Column-wise true/false value of a matrix
}

\description{
Column-wise true/false value of a matrix.
}

\usage{
colTrue(x)
colFalse(x)
colTrueFalse(x)
}

\arguments{
\item{x}{
A logical matrix with data.
}
}

%\details{
%
%}

\value{
An integer vector where item "i" is the number of the true/false values of "i" column.
}

%\references{
%Shabalin, A.A. (2012). Matrix eQTL: Ultra fast eQTL analysis %via large matrix operations. Bioinformatics 28(10): 1353-1358.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{rowMins}, \link{rowFalse}, \link{nth}, \link{colrange}, \link{colMedians}, \link{colVars}, \link{colSort}, \link{rowSort}, \link{rowTrue}
}
}

\examples{
x <- matrix(as.logical(rbinom(100*100,1,0.5)),100,100)

s1 <- colTrue(x) 

s1 <- colFalse(x)  

s1 <- colTrueFalse(x)

x<-s1<-NULL
}

\keyword{ Column-wise true }
\keyword{ column-wise false }
\keyword{ column-wise true-false }

