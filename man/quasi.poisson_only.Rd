\name{Many univariate simple quasi poisson regressions}
\alias{quasi.poisson_only}
\title{
Many univariate simple poisson regressions
}

\description{
It performs very many univariate simple poisson regressions.
}

\usage{
quasi.poisson_only(x, y, tol = 1e-09, maxiters = 100)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns are the variables. Currently only continuous variables are allowed.
}
\item{y}{
The dependent variable; a numerical vector with many discrete values (count data). 
}
\item{maxiters}{
The maximum number of iterations after which the Newton-Raphson algorithm is terminated.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
}

\details{
The function is written in C++ and this is why it is very fast. It can accept thousands of predictor variables. It is usefull for univariate screening. 
We provide no p-value correction (such as fdr or q-values); this is up to the user. 
}

\value{
A matrix with the deviance and the estimated phi parameter (dispersion parameter) of each simple poisson regression model for each predictor variable.
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Manos Papadakis <papadakm95@gmail.com> and Stefanos Fafalios <stefanosfafalios@gmail.com>

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>, Manos Papadakis <papadakm95@gmail.com> and Stefanos Fafalios <stefanosfafalios@gmail.com>.
}

\seealso{
\code{\link{poisson_only} \link{univglms}, \link{logistic_only}, \link{allbetas}, \link{regression}
}
}

\examples{
## 200 variables, hence 200 univariate regressions are to be fitted
x <- matrix( rnorm(100 * 200), ncol = 200 )
y <- rpois(100, 10)
system.time( poisson_only(x, y) )
b1 <- poisson_only(x, y) 
b2 <- quasi.poisson_only(x, y) 

b1<-b2<-x<-y<-NULL
}

\keyword{ Poisson regressions }
\keyword{ generalised linear models } 

