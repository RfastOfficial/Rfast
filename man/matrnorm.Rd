\name{Generates random values from a normal and puts them in a matrix}
\alias{matrnorm}
\title{
Generates random values from a normal and puts them in a matrix
}

\description{
Generates random values from a normal and puts them in a matrix.
}

\usage{
matrnorm(n, p, seed = NULL)
}

\arguments{
\item{n}{
The sample size, the number of rows the matrix will have. 
}
\item{p}{
The dimensionality of the data, the nubmer of columns of the matrix.
}
\item{seed}{
If you want the same to be generated again use a seed for the generator, an integer number.
}
}

\details{
How many times did you have to simulated data from a (standard) normal distribution in order to test 
something? For example, in order to see the speed of \code{\link{logistic_only}}, one needs to generate 
a matrix with predictor variables. The same is true for other similar functions. In \code{\link{sftests}},
one would like to examine the typer I error of this test under the null hypothesis. 

By using the Ziggurat method of generating standard normal variates, this function is really fast when you 
want to generate big matrices. 
}

\value{
An n x p matrix with data simulated from a standard normal distribution. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rvmf}, \link{Rnorm}, \link{rmvnorm}, \link{rvonmises}
}
}

\examples{
x <- matrnorm(100, 100)
}




