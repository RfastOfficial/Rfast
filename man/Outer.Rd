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
The available options are "*" (multiplication), "/" (division), "+" (sum), "-" (sbtraction), 
"^" (power raise), and "%%" (modulo operation). 
}
}

\details{
The function is the same as R's "outer", but works with vectors only and probably has less 
capabilities, but faster.
}

\value{
A matrix with two rows. In each row the X2 or G2 test statistic, its p-value and the degrees of freedom are returned.
}

%\references{
%
%}

\author{
Manos Papadakis and Michail Tsagris

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> 
and Michail Tsagris <mtsagris@yahoo.gr>.
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