\name{Gaussian regression with a log-link}
\alias{normlog.reg}
\title{
Gaussian regression with a log-link
}

\description{
Gaussian regression with a log-link.
}

\usage{
normlog.reg(y, x, tol = 1e-07, maxiters = 100)
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
A Gaussian regression with a log-link is fitted.
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
\item{be}{
The regression coefficients
}
}

\author{
Stefanos Fafalios

R implementation and documentation: Stefanos Fafalios <stefanosfafalios@gmail.com>
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{ \link{normlog.regs}, \link{score.glms}, \link{prop.regs}, \link{allbetas} }
}

\examples{
\donttest{
y <- abs( rnorm(100) )
x <- matrix( rnorm(100 * 2), ncol = 2)
a <- normlog.reg(y, x)
b <- glm(y ~ x, family = gaussian(log) )
summary(b)
a
}
}

