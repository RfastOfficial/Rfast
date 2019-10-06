\name{Equality of objects}
\alias{all_equals}
\title{Equality of objects
}

\description{Equality of objects.
}

\usage{
all_equals(x,y,round_digits = FALSE,without_attr=FALSE,fast_result=FALSE)
}

\arguments{
\item{x}{
A Matrix, List, Dataframe or Vector.
}
\item{y}{
A Matrix, List, Dataframe or Vector.
}
\item{round_digits}{
The digit for rounding numbers.
}
\item{without_attr}{
A boolean value (TRUE/FALSE) for deleting attributes. Be carefull although because some atributes are very important for you item.
}
\item{fast_result}{
A boolean value (TRUE/FALSE) for using just identical.But you can combine only with round_digits argument.
}
}

%\details{
%For every value of \eqn{\alpha} the normal likelihood (see the refernece) is computed. At the end, the plot of the values is constructed.
%}

\value{
A boolean (TRUE/FALSE) value which represents if the items x and y are equal.
}

%\references{
%Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{Match}, \link{mvbetas}, \link{correls}, \link{univglms}, \link{colsums}, \link{colVars}
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100 )
y <- matrix( rnorm(100 * 100), ncol = 100 )
all_equals(x,y) 
all_equals(x, x)
}

\keyword{ Equality check }

