\name{Some summary statistics of a vector for each level of a grouping variable}
\alias{group}
\alias{group.sum}
\alias{group.mean}
\title{
Some summary statistics of a vector for each level of a grouping variable. 
}

\description{
Some summary statistics of a vector for each level of a grouping variable.
}

\usage{
group(x,ina,method="sum",ina.min=NULL,ina.max = NULL,
	ina.length.unique=NULL,mad.method="median")
group.sum(x, ina,ina.max = NULL,ina.min = NULL)
group.mean(x, ina,ina.max = max(ina))
}

\arguments{
\item{x}{
A numerical vector with data.
}
\item{ina}{
A numerical vector with numbers. Note that zero and negative values are not allowed as this can cause R to run forever or crash.
}
\item{ina.length.unique}{
Length of the unique numerical values of ina argument.
}
\item{method}{
A character vector with values "sum", "var", "all", "any", "mad", "mean", "med", "min", "max", "min.max".
}
\item{ina.max}{
Maximum number for vector ina. 
}
\item{ina.min}{
Minimum number for vector ina. 
}
\item{mad.method}{
A character vector with values "median", for median absolute deviation or "mean", for mean absolute deviation. This works only with method="mad". 
}
}

\details{
This command works only for vectors. Median absolute deviation, mean, median, minimum, maximum are some of the options offered.
} 

\value{
A vector with the variance, or standard deviation, or mean, or minimum, or maximum, or median, or minimum-maximum of x for each distinct value of ina. 
}

\author{
Manos Papadakis and Michail Tsagris

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> and Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{colmeans}, \link{colVars}, \link{colMedians}
}
}

\examples{
\dontrun{
x <- rgamma(100,1, 4)
ina <- sample(1:5, 100, TRUE)
res<-group(x, ina,method="var")

}
}

\keyword{ Variance }
\keyword{ grouppings }
