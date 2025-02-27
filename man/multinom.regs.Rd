\name{Many simple multinomial regressions}
\alias{multinom.regs}
\title{
Many simple multinomial regressions
}

\description{
Many simple multinomial regressions.
}

\usage{
multinom.regs(y, x, tol = 1e-08, logged = FALSE, parallel = FALSE, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, either a numerical variable or a factor variable.
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
Do you want this to be executed in parallel or not. The parallel takes place in C++, and the number of threads 
is defined by each system's availiable cores.
}
\item{maxiters}{
The maximum number of iterations that can take place in each regression.
}
}

\details{
Many simple multinomial regressions are fitted.
}

\value{
A matrix with the test statistic values, their relevant (logged) p-values and the BIC values.
}

\author{
Stefanos Fafalios

R implementation and documentation: Stefanos Fafalios <stefanosfafalios@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{poisson_only}, \link{prop.regs}, \link{score.geomregs} }
}

\examples{
y <- rbinom(100, 2, 0.5)
x <- matrnorm(100, 100)
a <- multinom.regs(y, x)
x <- NULL
}

