\name{Spatial median regression}
\alias{spatmed.reg}
\title{
Spatial median regression
}

\description{
Spatial median regression with Euclidean data.
}

\usage{
spatmed.reg(y, x, tol = 1e-07)
}

\arguments{
\item{y}{
A matrix with the response variable. 
}
\item{x}{
The predictor variable(s), they have to be continuous.
}
\item{tol}{
The threshold upon which to stop the iterations of the Newton-Rapshon algorithm.
}
}

\details{
The objective function is the minimization of the sum of the absolute residuals. 
It is the multivariate generalisation of the median regression.
}

\value{
A list including:
\item{iters}{
The number of iterations that were required.
}
\item{be}{
The beta coefficients.
}
}

\references{
Biman Chakraborty (2003) On multivariate quantile regression. Journal of Statistical 
Planning and Inference http://www.stat.nus.edu.sg/export/sites/dsap/research/documents/tr01_2000.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{spat.med}, \link{sscov}, \link{lmfit}
}
}

\examples{
\donttest{
x <- as.matrix(iris[, 3:4])
y <- as.matrix(iris[, 1:2])
mod1 <- spatmed.reg(y, x)
}
}


