\name{Many simple regressions for positive valued data}
\alias{normlog.regs}
\alias{gammaregs}
\alias{invgauss.regs}
\title{
Many simple regressions for positive valued data
}

\description{
Many simple regressions for positive valued data.
}

\usage{
normlog.regs(y, x, tol = 1e-08, logged = FALSE, parallel = FALSE, maxiters = 100)
gammaregs(y, x, tol = 1e-07, logged = FALSE, maxiters = 100)
invgauss.regs(y, x, tol = 1e-08, logged = FALSE, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical variable with non negative numbers for the Gamma and inverse Gaussian regressions. 
For the Gaussian with a log-link zero values are allowed.
}
\item{x}{
A matrix with the indendent variables.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
\item{parallel}{
Do you want this to be executed in parallel or not. The parallel takes place in C++, therefore you do not have the option to set the number of cores.
}
\item{maxiters}{
The maximum number of iterations that can take place in each regression.
}
}

\details{
Many simple Gamma, inverse Gaussian or Gaussian regressions with a log-link are fitted.
}

\value{
A matrix with the test statistic values and their relevant (logged) p-values.
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989. 

Zakariya Yahya Algamal and Intisar Ibrahim Allyas (2017). Prediction of blood lead level in maternal and fetal 
using generalized linear model. International Journal of Advanced Statistics and Probability, 5(2): 65-69.
}

\author{
Stefanos Fafalios and and Michail Tsagris

R implementation and documentation: Stefanos Fafalios <stefanosfafalios@gmail.com> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{normlog.reg}, \link{score.glms}, \link{prop.regs}, \link{allbetas} }
}

\examples{
y <- abs( rnorm(100) )
x <- matrnorm(100, 100)
a <- normlog.regs(y, x)
b <- glm(y ~ x[, 1], family = gaussian(log) )
anova(b, test= "F")
a[1, ]
a2 <- gammaregs(y, x)
a3 <- invgauss.regs(y, x)
x <- NULL
}

