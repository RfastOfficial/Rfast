\name{Convert a dataframe to matrix}
\alias{data.frame.to_matrix}
\title{
Convert a dataframe to matrix
}

\description{
Convert a dataframe to matrix.
}

\usage{
data.frame.to_matrix(x,col.names = NULL,row.names = NULL)
}

\arguments{
\item{x}{
A Numeric matrix with data and NAs.
}
\item{col.names}{
A boolean value for keeping the colnames for argument x or a character vector for the new colnames.
}
\item{row.names}{
A boolean value for keeping the rownames for argument x or a character vector for the new rownames.
}
}

\details{
This functions converts a dataframe to matrix. Even if there are factors, the function converts them into numerical values. Attributes are not allowed for now.
}

\value{
A matrix wich has the numrical values from the dataframe.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>
}


\seealso{
\code{\link{Match}, \link{is.symmetric}, \link{permutation}
}
}

\examples{
res<-data.frame.to_matrix(iris)
}

\keyword{ Dataframe to Matrix }