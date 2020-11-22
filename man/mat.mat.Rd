\name{Number of equal columns between two matrices}
\alias{mat.mat}
\title{Number of equal columns between two matrices
}

\description{Number of equal columns between two matrices.
}

\usage{
mat.mat(x, y)
}

\arguments{
\item{x}{
A numerical matrix. See details for more information. It must have the same number of rows as y.
}
\item{y}{
A numerical matrix. See details for more information. It must have the same number of rows as x.
}
}

\details{
The function takes each column of x and checks the number of times it matches a column of y. In the example below, we take the first 3 columns of iris 
as the x matrix. The y matrix is the whole of iris. We will see how many times, each column of x appears in the y matrix. The answer is 1 for each column. 

}

\value{
A numerical vector of size equal to the number of columns of x. 
}

\author{
Manos Papadakis 

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{Match}, \link{colmeans}, \link{colMedians}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
y <- iris
y[, 5] <- as.numeric(y[, 5])
y <- as.matrix(y)
res<-mat.mat(x, y)

x<-y<-NULL
}

\keyword{ Matrices }

