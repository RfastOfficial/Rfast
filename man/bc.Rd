\name{Estimation of the Box-Cox transformation}
\alias{bc}
\title{
Estimation of the Box-Cox transformation
}

\description{
Estimation of the Box-Cox transformation.
}

\usage{
bc(x, low = -1, up = 1)
}

\arguments{
\item{x}{
A numerical vector with strictly positive values.
}
\item{low}{
The lowest value to search for the best \eqn{\lambda} parameter.
}
\item{up}{
The highest value to search for the best \eqn{\lambda} parameter.
}
}

\details{
The functions estimates the best \eqn{\lambda} in the Box-Cox power transformation.
}

\value{
The optimal value of \eqn{\lambda}.
}

\references{
Box George E. P. and Cox D. R. (1964). An analysis of transformations. 
Journal of the Royal Statistical Society, Series B, 26 (2):211-252. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{correls}, \link{auc}
}
}

\examples{
x <- exp(rnorm(1000))
res<-bc(x)
}


