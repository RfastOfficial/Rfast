\name{Correlations}
\alias{correls}
\alias{groupcorrels}
\title{
Correlation between a vector and a set of variables
}

\description{
Correlation between a vector and a set of variables.
}

\usage{
correls(y, x, type = "pearson", a = 0.05, rho = 0)
groupcorrels(y, x, type = "pearson", ina)
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A matrix with the data.
}
\item{type}{
The type of correlation you want. "pearson" and "spearman" are the two supported types for the "correls" because their standard error is easily calculated. For the "groupcorrels" you can also put "kendall" because no hypothesis test is performed in that function.
}
\item{a}{
The significance level used for the confidence intervals.
}
\item{rho}{
The value of the hypothesised correlation to be used in the hypothesis testing.
}
\item{ina}{
A factor variable or a numeric variable idicating the group of each observation.
}
}

\details{
The functions uses the built-in function "cor" which is very fast and then includes confidence intervals and produces a p-value for the hypothesis test. 
}

\value{
For the "correls" a matrix with 5 column; the correlation, the p-value for the hypothesis test that each of them is eaqual to "rho", the test statistic and the $a/2\%$ lower and upper confidence limits.

For the "groupcorrels" a matrix with rows equal to the number of groups and columns equal to the number of columns of x. The matrix contains the correlations only, no statistical hypothesis test is performed.
}

%\references{
%Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{allbetas}, \link{univglms}
}
}

\examples{
x <- matrnorm(60, 100 )
y <- rnorm(60)
r <- cor(y, x)  ## correlation of y with each of the xs
a <- allbetas(y, x)  ## the coefficients of each simple linear regression of y with x
b <- correls(y, x)
ina <- rep(1:2, each = 30)
b2 <- groupcorrels(y, x, ina = ina) 
x <- NULL
}

\keyword{ Correlations }


