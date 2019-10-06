\name{Row-wise true value }
\alias{rowTrue}
\alias{rowFalse}
\alias{rowTrueFalse}
\title{
Row-wise true value of a matrix
}

\description{
Row-wise true value of a matrix.
}

\usage{
rowTrue(x)
rowFalse(x)
rowTrueFalse(x)
}

\arguments{
\item{x}{
A logical matrix with data.
}
}


\value{
An integer vector where item "i" is the number of the true/false values of "i" row.
}

%\references{
%Shabalin, A.A. (2012). Matrix eQTL: Ultra fast eQTL analysis %via large matrix operations. Bioinformatics 28(10): 1353-1358.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{rowMins}, \link{colFalse}, \link{nth}, \link{rowrange}, \link{rowMedians}, \link{rowVars}, \link{colTrue}
}
}

\examples{
x <- matrix(as.logical(rbinom(100*100,1,0.5)),100,100)

s1 <- rowTrue(x) 

s1 <- rowFalse(x)  

s1 <- rowTrueFalse(x)

x<-s1<-NULL
}

\keyword{ Row-wise true }
\keyword{ Row-wise false }
\keyword{ Row-wise true-false }

