\name{Many univariate simple logistic and Poisson regressions}
\alias{logistic_only}
\alias{poisson_only}

\title{
Many univariate simple binary logistic regressions
}

\description{
It performs very many univariate simple binary logistic regressions.
}

\usage{
logistic_only(x, y, tol = 1e-09, b_values = FALSE)
poisson_only(x, y, tol = 1e-09, b_values = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns are the variables. Currently only continuous variables are allowed.
}
\item{y}{
The dependent variable; a numerical vector with two values (0 and 1) for the logistic regressions and 
a vector with many discrete values (count data) for the Poisson regressions. 

}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{b_values}{
Do you want the values of the coefficients returned? If yes, set this to TRUE.
}
}

\details{
The function is written in C++ and this is why it is very fast. It can accept thousands of predictor variables. It is usefull for univariate screening. 
We provide no p-value correction (such as fdr or q-values); this is up to the user. 
}

\value{
A vector with the deviance of each simple binayr logistic regression model for each predictor variable.
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Manos Papadakis <papadakm95@gmail.com>

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{univglms}, \link{score.glms}, \link{prop.regs}, \link{quasi.poisson_only}, \link{allbetas}, \link{correls}, \link{regression}
}
}

\examples{
\dontrun{
## 300 variables, hence 300 univariate regressions are to be fitted
x <- matrix( rnorm(100 * 300), ncol = 300 )

## 100 observations in total
y <- rbinom(100, 1, 0.6)   ## binary logistic regression
system.time( logistic_only(x, y) )
a1 <- logistic_only(x, y)
 
a2 <- glm(y ~ x[, 1], binomial)$deviance 
a2 <- as.vector(a2)

y <- rpois(100, 10)
system.time( poisson_only(x, y) )

a1 <- x <- NULL
}
}

\keyword{ Logistic regressions }
\keyword{ generalised linear models } 

