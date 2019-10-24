\name{Inverse Gaussian regression with a log-link}
\alias{invgauss.reg}
\title{
Inverese Gaussian regression with a log-link
}

\description{
Inverse Gaussian regression with a log-link.
}

\usage{
invgauss.reg(y, x, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical variable with non negative numbers.
}
\item{x}{
A matrix or data.frame with the indendent variables.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations that can take place in the regression.
}
}

\details{
An inverse Gaussian regression with a log-link is fitted.
}

\value{
A list including:
\item{i}{
The number of iterations required by the Newton-Raphson
}
\item{loglik}{
The log-likelihood value.
}
\item{deviance}{
The deviance value.
}
\item{phi}{
The dispersion parameter (\eqn{\phi}) of the regression. This is necessary if you want to perform an F hypothesis test for 
the significance of one or more independent variables.
}
\item{be}{
The regression coefficients
}
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989. 

Zakariya Yahya Algamal and Intisar Ibrahim Allyas (2017). Prediction of blood lead level in maternal and fetal 
using generalized linear model. International Journal of Advanced Statistics and Probability, 5(2): 65-69.
}

\author{
Michail Tsagris

R implementation and documentation: Stefanos Fafalios <mtsagris@yahoo.gr>
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{ \link{invgauss.regs}, \link{normlog.reg}, \link{score.glms} }
}

\examples{
\dontrun{
y <- abs( rnorm(100) )
x <- matrix( rnorm(100 * 2), ncol = 2)
a <- invgauss.reg(y, x)
a
}
}

