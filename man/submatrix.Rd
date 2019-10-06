\name{Sub-matrix}
\alias{submatrix}
\title{Sub-matrix
}

\description{Sub-matrix.
}

\usage{
submatrix(x,rowStart=1,rowEnd=1,colStart=1,colEnd=1)
}

\arguments{
\item{x}{
A Matrix, List, Dataframe or Vector.
}
\item{rowStart}{
Start of the row.
}
\item{rowEnd}{
End of the row.
}
\item{colStart}{
Start of the col.
}
\item{colEnd}{
End of the col.
}
}

\value{
sub matrix like R's, x[startrow:endrow,startcol:endcol]. Fast especially for big sub matrices.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{ \link{Match}, \link{mvbetas}, \link{correls}, \link{univglms}, \link{colsums}, \link{colVars}
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100 )
res<-submatrix(x,1,50,1,25) # x[1:50,1:25]

x<-NULL
}

\keyword{ Sub-matrix }

