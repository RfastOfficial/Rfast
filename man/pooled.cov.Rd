\name{Pooled covariance matrix}
\alias{pooled.cov}
\title{
Pooled covariance matrix
}

\description{
Pooled covariance matrix.
}

\usage{
pooled.cov(x, ina)
}

\arguments{
\item{x}{
A matrix with continuous data.
}
\item{ina}{
A numerical vector indicating the groups. \bold{The nubmers must be consecutive and start from 1}.
}
}

\details{
The spatial median is at first computed (if not supplied) and then the covariance matrix. 
}

\value{
The spatial sign covariance matrix.
}

\references{
Durre A, Vogel D. and D.E. Tyler D.E.(2014). The spatial sign covariance matrix with unknown location.  Journal of Multivariate Analysis, 130: 107-117.
http://arxiv.org/pdf/1307.5706v2.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spat.med}, \link{spatmed.reg}
}
}

\examples{

res<-sscov( as.matrix(iris[, 1:4]) )

}

\keyword{ covariance matrix }
\keyword{ robust statistics }
