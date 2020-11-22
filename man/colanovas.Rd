\name{Many Welch's F-tests}
\alias{colanovas}
\title{
Many Welch's F-tests
}

\description{
Many Welch's F-tests.
}

\usage{
colanovas(y, x, logged = FALSE)
}

\arguments{
\item{y}{
A numerical vector with the dependent variable. 
}
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the 
columns are the variables. This must be a matrix with the categorical variables as numbers, starting from 1.
Welch's F-test is performed for each variable. 
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
}

\details{
For each categorical variable in the x matrix Welch's F test is performed. This is the opposie of \code{\link{ftests},
where there are many dependent variables and one categorical variable.
}
}

\value{
A matrix with the test statistic and the p-value for each predictor variable.
}

\references{
Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{regression}, \link{ftests}, \link{allbetas}, \link{correls}
}
}

\examples{
y <-  rnorm(100)
x <- matrix( rbinom(100 * 50, 2, 0.5) + 1 , ncol = 50)  
a <- colanovas(y, x)
x <- NULL
}



