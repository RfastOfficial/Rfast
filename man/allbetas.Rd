\name{Many simple linear regressions coefficients}
\alias{allbetas}
\title{Simple linear regressions coefficients
}

\description{Simple linear regressions coefficients.
}

\usage{
allbetas(y, x, pvalue = FALSE, logged = FALSE)
}

\arguments{
\item{y}{
A numerical vector with the response variable. 
}
\item{x}{
A matrix with the data, where rows denotes the observations and the columns contain the independent variables.
}
\item{pvalue}{
If you want a hypothesis test that each slope (beta coefficient) is equal to zero set this equal to TRUE. It will also produce all the correlations between y and x.
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
}

%\details{
%For every value of \eqn{\alpha} the normal likelihood (see the refernece) is computed. At the end, the plot of the values is constructed.
%}

\value{
A matrix with the constant (alpha) and the slope (beta) for each simple linear regression. 
If the p-value is set to TRUE, the correlation of each y with the x is calculated along with the relevant test statistic and its associated p-value.
}

%\references{
%Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{mvbetas}, \link{correls}, \link{univglms}, \link{colsums}, \link{colVars}
}
}

\examples{
x <- matrix( rnorm(100 * 50), ncol = 50 )
y <- rnorm(100)
r <- cor(y, x)  ## correlation of y with each of the xs
a <- allbetas(y, x)  ## the coefficients of each simple linear regression of y with x
x <- NULL
}

\keyword{ Many betas in regression }

