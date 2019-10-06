\name{Many multivariate simple linear regressions coefficients}
\alias{mvbetas}
\title{Many multivariate simple linear regressions coefficients
}

\description{Many multivariate simple linear regressions coefficients.
}

\usage{
mvbetas(y, x, pvalue = FALSE)
}

\arguments{
\item{y}{
A matrix with the data, where rows denotes the observations and the columns contain the dependent variables.
}
\item{x}{
A numerical vector with one continuous independent variable only. 
}
\item{pvalue}{
If you want a hypothesis test that each slope (beta coefficient) is equal to zero set this equal to TRUE. It will also produce all the correlations between y and x.
}
}

\details{
It is a function somehow opposite to the \code{\link{allbetas}}. Instead of having one y and many xs we have many ys and one x.
}

\value{
A matrix with the constant (alpha) and the slope (beta) for each simple linear regression. 
If the p-value is set to TRUE, the correlation of each y with the x is calculated along with the relevant p-value.
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{allbetas}, \link{correls}, \link{univglms}
}
}

\examples{
y <- matrnorm(100, 100)
x <- rnorm(100)
a <- mvbetas(y, x, pvalue = FALSE)
b <- matrix(nrow = 100, ncol = 2)
z <- cbind(1, x)

system.time( a <- mvbetas(y, x) )
system.time(  for (i in 1:100) b[i, ] = coef( lm.fit( z, y[, i] ) )  )
y <- NULL
}

\keyword{ Many betas in regression }

