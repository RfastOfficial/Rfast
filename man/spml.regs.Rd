\name{Many simple circular or angular regressions}
\alias{spml.regs}
\title{
Many simple circular or angular regressions
}

\description{
Many regressions with one circular dependent variable and one Euclidean independent variable.
}

\usage{
spml.regs(y, x, tol = 1e-07, logged = FALSE, maxiters = 100, parallel = FALSE)
}

\arguments{
\item{y}{
The dependent variable, it can be a numerical vector with data expressed in radians or it can be a matrix with
two columns, the cosinus and the sinus of the circular data. The benefit of the matrix is that if the function
is to be called multiple times with the same response, there is no need to transform the vector every time into
a matrix.
}
\item{x}{
A matrix with independent variable.
}
\item{tol}{
The tolerance value to terminatate the Newton-Raphson algorithm.
}
\item{logged}{
Do you want the logarithm of the p-value be returned? TRUE or FALSE.
}
\item{maxiters}{
The maximum number of iterations to implement.
}
\item{parallel}{
Do you want the calculations to take plac ein parallel? The default value if FALSE.
}
}

\details{
The Newton-Raphson algorithm is fitted in these regression as described in Presnell et al. (1998). For each colum of x a circual
regression model is fitted and the hypothesis testing of no association between y and this variable is performed.
}

\value{
A matrix with two columns, the test statistics and their associated (log) p-values.
}

\references{
Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for
directional data. Journal of the American Statistical Association, 93(443): 1068-1077.
}

\author{
Michail Tsagris and Stefanos Fafalios

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Stefanos Fafalios
<stefanosfafalios@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spml.mle}, \link{iag.mle}, \link{acg.mle}
}
}

\examples{
x <- rnorm(100)
z <- cbind(3 + 2 * x, 1 -3 * x)
y <- cbind( rnorm(100,z[ ,1], 1), rnorm(100, z[ ,2], 1) )
y <- y / sqrt( rowsums(y^2) )
x <- matrnorm(100, 100)
a <- spml.regs(y, x)
x <- NULL
}

\keyword{ Circular regression }
\keyword{ projected normal }
