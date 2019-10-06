\name{Many simple geometric regressions}
\alias{geom.regs}
\title{
Many simple geometric regressions.
}

\description{
Many simple geometric regressions.
}

\usage{
geom.regs(y, x, tol = 1e-07, type = 1, logged = FALSE, parallel = FALSE, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, count data.
}
\item{x}{
A matrix with the indendent variables.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{type}{
Type 1 refers to the case where the minimum is zero and type 2 for the case of the minimum being 1.
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
\item{parallel}{
Do you want this to be executed in parallel or not. The parallel takes place in C++, and the number of threads is defined by each system's availiable cores.
}
\item{maxiters}{
The max number of iterations that can take place in each regression.
}
}

\details{
Many simple geometric regressions are fitted.
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
y <- rgeom(100, 0.6)
x <- matrix( rnorm(100 * 50), ncol = 50)
a <- geom.regs(y, x)
x <- NULL
}

