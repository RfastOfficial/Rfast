\name{Column and row wise coefficients of variation}
\alias{colcvs}
\alias{rowcvs}
\title{
Column and row wise coefficients of variation
}

\description{
Column and row wise coefficients of variation.
}

\usage{
colcvs(x, ln = FALSE, unbiased = FALSE) 
rowcvs(x, ln = FALSE, unbiased = FALSE) 
}

\arguments{
\item{x}{
A numerical matrix with the data. 
}
\item{ln}{
If you have log-normally distributed data (or assume you do), then set this to TRUE. 
}
\item{unbiased}{
A boolean variable indicating whether the unbiased for shpould be returned. This is applicable in case of small samples.
}
}

\details{
The colum-wise coefficients of variation are calculated.
}

\value{
A vector with the coefficient of variation for each column or row.
}

%\references{
%Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons 
%(excellent book).
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{colsums}, \link{colVars}
}
}

\examples{
m <- rnorm(100, 10)
x <- matrix(rnorm(100 * 100, m, 1), ncol = 100)
a1 <- colcvs(x)
a2 <- colcvs(x[1:25, ], unbiased = TRUE)
a3 <- colcvs( exp(x), ln = TRUE)
x <- NULL
}

\keyword{ Linear mixed models }
\keyword{ balanced design }
\keyword{ maximum likelihood estimation }
