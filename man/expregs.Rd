\name{Many exponential regressions}
\alias{expregs}
\title{
Many exponential regressions}

\description{
Many exponential regressions.
}

\usage{
expregs(y, x, di, tol = 1e-09, logged = FALSE)  
}

\arguments{
\item{y}{
A vector with positive data (including zeros).
}
\item{x}{
A numerical matrix with the predictor variables.  
}
\item{di}{
A vector of size equal to that of y with 0s and 1s indicating censoring or not respectively.
}
\item{tol}{
The tolerance value to stop the newton-Raphson iterations. It is set to 1e-09 by default. 
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
}

\details{
We have implemented the newton-Raphson in order to avoid unnecessary calculations.  
}

\value{
A matrix with three columns, the test statistic, its associated (logged) p-value and the BIC of each model.   
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{univglms}, \link{score.glms}, \link{logistic_only},  \link{poisson_only}, \link{regression}
}
}

\examples{
## 200 variables, hence 200 univariate regressions are to be fitted
x <- matrnorm(100, 100)
y <- rexp(100, 4)
system.time( expregs(y, x, di = rep(1, length(y))) )
x <- NULL
}

\keyword{ Exponential regressions }
