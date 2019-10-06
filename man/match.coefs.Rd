\name{Column-wise matching coefficients}
\alias{match.coefs}
\title{
Column-wise matching coefficients
}

\description{
Column-wise matching coefficients.
}

\usage{
match.coefs(x, y = NULL, ina, type = "jacc") 
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples and the columns are the variables.
}
\item{y}{
A second matrix with the data of the second group. If this is NULL (default value) then the argument ina must be 
supplied. Notice that when you supply the two matrices the procedure is two times faster.
} 
\item{ina}{
A numerical vector with 1s and 2s indicating the two groups. Be careful, the function is designed to accept only
these two numbers. In addition, if your "y" is NULL, you must specify "ina".
}
\item{type}{
This denotes the type of matching coefficient to calculate. For the Jaccard index put "jacc". For the simple 
matching coefficient put "smc" or else both of them will be calculated. 
}
}

\details{
Two matrices are given as imput and for each column matching coefficients are calculated, either the Jaccard or 
the simple matching coefficient or both.  
}

\value{
A matrix with one or two columns, depending on the type you have specified. If you specify "both", there will be
two columns, if you specify "jacc" or "smc" then just one column.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{odds}, \link{colTabulate}
}
}

\examples{
x <- matrix(rbinom(400 * 10, 1, 0.5), ncol = 10)
y <- matrix(rbinom(400 * 10, 1, 0.5), ncol = 10)
a <- match.coefs(x, y, type = "both")
x <- NULL
y <- NULL
}


