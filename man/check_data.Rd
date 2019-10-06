\name{Search for variables with zero range in a matrix}
\alias{check_data}
\title{Search for variables with zero range in a matrix
}

\description{Search for variables with zero range in a matrix.
}

\usage{
check_data(x, ina = NULL) 
}

\arguments{
\item{x}{
A matrix or a data.frame with the data, where rows denotes the observations and the columns contain the dependent variables.
}
\item{ina}{
If your data are grouped, for example there is a factor or numerical variable indicating the groups of the data supply it here, otherwise leave it NULL. 
}
}

\details{
The function identifies the variabels with zero range, instead of a zero variance as this is faster. It will work with matrices and data.frames. 
}

\value{
A numerical vector of length zero if no zero ranged variable exists, or of length at least one with the index (or indices) of the variable(s) that need 
attention or need to be removed.
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{colrange}, \link{colVars}
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100 )
res<-check_data(x)

## some variables have a constant value
x[, c(1,10, 50, 70)] <- 1
res<-check_data(x)
id <- rep(1:4, each = 25 )
x[1:25, 2] <- 0
res<-check_data(x)  ## did not use the id variable
res<-check_data(x, id)  ## see now
x <- NULL
}

\keyword{ Zero range }
\keyword{ data check }
