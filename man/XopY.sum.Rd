\name{Operations between two matrices or matrix and vector}
\alias{XopY.sum}
\alias{eachrow}
\alias{eachcol.apply}
\title{
Operations between two matrices or matrix and vector
}

\description{
Operations between two matrices or matrix and vector.
}

\usage{
XopY.sum(x, y = NULL, oper = "*")
eachrow(x,y,oper = "*",method = NULL)
eachcol.apply(x,y,indices = NULL,oper = "*",apply = "sum", parallel = FALSE)
}

\arguments{
\item{x}{
A numerical matrix.
}
\item{y}{
A second numerical matrix for "XopY.sum" whose dimensions must match the ones of x, or vector for "eachrow","eachcol.apply" whose length must match with the rows of x.
}
\item{oper}{
The operation to be performed, either "*", "/", "+", "-" or "==".
}
\item{method}{
A character value for choosing option to apply in the result.
Options:
	1) sum
	2) max
	3) min

Does not work for oper="==".

}
\item{indices}{
An integer vector with indices to specific columns. Only for "eachcol.apply".
}
\item{apply}{
A character value with the function to be applied in the columns of the matrix. Only for "eachcol.apply". 
Options:
	1) sum
	2) median
	3) max
	4) min
}
\item{parallel}{
A boolean value for parallel version.
}

\details{
XopY.sum: sum(X op Y) where op can be on of "+,-,*,/".

eachrow: X op Y \bold{by row} or \bold{FUNCTION}(X op Y) where "x" is matrix, "y" is vector with length as much an the columns of x and "op" is one of "+,-,*,/,==", and "FUNCTION" is a specific method for applying in the result matrix (see argument method).

eachcol.apply: \bold{FUNCTION}(X op Y) \bold{by column} where "x" is matrix, "y" is vector with length as much an the rows of x, "op" is one of "+,-,*,/" and "FUNCTION" is a specific method (see argument apply).

\bold{NOTE:} Arguments "method" does not work for oper="==" and this operation works only in "eachrow".
}

\value{
XopY.sum: sum(X op Y) where "op" can be on of "+,-,*,/".

eachrow: operation by row between a matrix and a vector."op" can be on of "+,-,*,/".
If "suma=TRUE" then returns the sum of this operation.

eachcol.apply: operation by column between a matrix and a vector and applied a specific function."op" can be on of "+,-,*,/".
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{Dist}, \link{dista}, \link{colmeans}, \link{Diag.fill}, 
\link{colMads}, \link{rowMads}
}
}

\examples{
x <- matrix( rnorm(5 * 5), ncol = 5 )
y <- matrix( rnorm(5 * 5), ncol = 5 )
res<-XopY.sum(x, y, oper = "*")
y <- x[,1]
res<-eachrow(x,y)

all.equal(eachcol.apply(x,y),colsums(x*y))


x<-y<-NULL
}

\keyword{ Sum }
\keyword{ operations }