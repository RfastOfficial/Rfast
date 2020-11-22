\name{Column-wise Yule's Y (coefficient of colligation)}
\alias{col.yule}
\title{
Column-wise Yule's Y (coefficient of colligation)
}

\description{
Column-wise Yule's Y (coefficient of colligation).
}

\usage{
col.yule(x, y = NULL, ina)
}

\arguments{
\item{x}{
A matrix with 0 and 1. Every column refers to a different sample or variable. 
}
\item{y}{
A second matrix, of the same dimensions as x, with 0 and 1. Every column refers to a different sample or variable. 
}
\item{ina}{
If y is NULL, ina must be specified. This is a numeric vector with 1s and 2s, indicating the group of each row.
}
}

\details{
Yule's coefficient of colligation is calculated for every column.
}

\value{
A vector with Yule's Y, one for every column of x is returned.
}

\references{
Yule G. Udny (1912). On the Methods of Measuring Association Between Two Attributes. 
Journal of the Royal Statistical Society, 75(6):579-652.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{yule}, \link{odds}
}
}

\examples{
x <- matrix(rbinom(300 * 10, 1, 0.5), ncol = 10)
ina <- rep(1:2, each = 150)
res<-col.yule( x, ina = ina )
}


