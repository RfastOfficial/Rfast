\name{Many odds ratio tests}
\alias{odds}
\title{
Many odds ratio tests
}

\description{
It performs very many odds ratio tests.
}

\usage{
odds(x, y = NULL, ina, logged = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the observations and the columns are the variables. 
They must be 0s and 1s only.
}
\item{y}{
A second matrix with the data of the second group. If this is NULL (default value) then the argument ina must 
be supplied. Notice that when you supply the two matrices the procedure is two times faster. They must be 0s and 1s
only.
} 
\item{ina}{
A numerical vector with 1s and 2s indicating the two groups. Be careful, the function is designed to accept only these two numbers. 
In addition, if your "y" is NULL, you must specify "ina".
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
Many odds ratio tests are performed. 
}

\value{
A matrix with the test statistic and the p-value (or their logarithm) of each test.
}

\references{
Mosteller Frederick (1968). Association and Estimation in Contingency Tables. Journal of the American 
Statistical Association. 63(321):1-28.
 
Edwards A.W.F. (1963). The measure of association in a 2x2 table. Journal of the Royal Statistical Society, Series A. 
126(1):109-114. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{odds.ratio}, \link{g2Test_univariate}
}
}

\examples{
x <- matrix( rbinom(100 * 100, 1, 0.5), ncol = 100 )
ina <- rep(1:2, each = 50)
a <- odds(x, ina = ina)
}

\keyword{ Odds ratios }


