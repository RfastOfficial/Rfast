\name{Outer function}
\alias{Outer}
\title{
Outer function
}

\description{
The outer function.
}

\usage{
Outer(x, y, oper = "*")
}

\arguments{
\item{x}{
A numerical vector.
}
\item{y}{
A numerical vector.
}
\item{oper}{
The available options are "*" (multiplication), "/" (division), "+" (sum), "-" (substraction), 
"^" (power raise), and "%%" (modulo operation). 
}
}

\details{
The function is the same as R's "outer", but works with vectors only and probably has less 
capabilities, but faster.
}

\value{
A matrix with all the combinations.
}

%\references{
%
%}

\author{
Manos Papadakis and Michail Tsagris

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> 
and Michail Tsagris <mtsagris@uoc.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{mat.mult}, \link{vecdist}
}
}

\examples{
x <- rnorm(10)
y <- rnorm(10)
res<-Outer(x, y)
}
