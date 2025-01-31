\name{Distance variance, covariance and correlation}
\alias{dvar}
\alias{dcov}
\alias{dcor}
\alias{bcdcor}
\title{
Distance variance, covariance and correlation
}

\description{
Distance variance, covariance and correlation.
}

\usage{
dvar(x, bc = FALSE)
dcov(x, y, bc = FALSE)
dcor(x, y, bc = FALSE)
bcdcor(x, y)
}

\arguments{
\item{x}{
A numerical matrix or a vector.
}
\item{y}{
A numerical matrix or a vector.
}
\item{bc}{
If you want the bias-corrected distance correlation set this equal to TRUE.
}
}

\details{
The distance variance of a matrix/vector, the distance covariance or distance correlation of two matrices is calculated.
For the distance variance of a vector we use the fast method of Huo and Szekely (2016).
}

\value{
The distance covariance or the distance variance.

For the bias corrected distance correlation its value only.
 
For the distance correlation a list including:
\item{dcov}{
The distance covariance.
}
\item{dvarX}{
The distance variance of x.
}
\item{dvarY}{
The distance variance of Y.
}
\item{dcor}{
The distance correlation.
}
}

\references{
Szekely G.J., Rizzo M.L. and Bakirov N.K.(2007). Measuring and Testing Independence 
by Correlation of Distances. Annals of Statistics, 35(6):2769--2794.

Szekely G. J. and Rizzo M. L. (2023). The Energy of Data and Distance Correlation. 
Chapman and Hall/CRC.

Tsagris M. and Papadakis M. (2025). Fast and light-weight energy statistics using the R package Rfast.
https://arxiv.org/abs/2501.02849
}

\author{
Michail Tsagris and Manos Papadakis.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{pdcor}, \link{dcor.ttest}, \link{edist}
}
}

\examples{
x <- as.matrix(iris[1:50, 1:4])
y <- as.matrix(iris[51:100, 1:4])
res <- dvar(x[, 1])
#dcor(x, y)
}

\keyword{ Distance variance }
\keyword{ Distance covariance }
