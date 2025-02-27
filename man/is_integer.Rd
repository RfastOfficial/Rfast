\name{Check if values are integers and convert to integer}
\alias{is_integer}
\alias{as_integer}
\title{
Check if values are integers and convert to integer
}

\description{
Check if values are integers and convert to integer.
}

\usage{
is_integer(x)
as_integer(x,result.sort = TRUE,init = 1, parallel = FALSE)
}

\arguments{
\item{x}{
is_integer: A vector with numeric data.
as_integer: A vector with data.
}
\item{result.sort}{
A logical value for sorting the result.
}
\item{init}{
An integer value to start.
}
\item{parallel}{
Run the algorithm in parallel. 
}
}

\details{
The behavior of these functions are different than R's built in. 

is_integer: 
	check if all the values are integers in memory. If typeof is double, and the values are integers in range -2^31 : 2^31 then it is  better to convert to integer vector for using less memory. Also you can decrease the time complexity.

as_integer: 
	converts the discrete values to integers.

}

\value{
is_integer: A logical value, TRUE if all values are integers and in range -2^31 : 2^31. Otherwise FALSE.

as_integer: By default the function will return the same result with "as.numeric" but the user can change the "init" value not start from 1 like R's. Also the result can be unsorted using "result.sort".
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{as_integer}, \link{colVars}, \link{colmeans}
}
}

\examples{
x<-runif(10)
y1<-is_integer(x) # y1 is FALSE
x<-as.numeric(rpois(10,10)) # integers but typeof is double
y1<-is_integer(x) # y1 is TRUE so you can convert to integer vector.

as_integer(letters) ## as.numeric(letters) produce errors
x<-y1<-NULL
}

\keyword{ is_integer Creation }