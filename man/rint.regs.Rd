\name{Many simple linear mixed model regressions}
\alias{rint.regs}
\title{
Many simple linear mixed model regressions
}

\description{
Many simple linear mixed model regressions with random intercepts only.
}

\usage{
rint.regs(y, x, id, tol = 1e-08, logged = FALSE, parallel = FALSE, maxiters = 100) 
}

\arguments{
\item{y}{
A numerical vector with the data. The subject values, the clustered data. 
}
\item{x}{
A numerical matrix with data ,the independent variables. 
}
\item{id}{
A numerical variable with 1, 2, ... indicating the subject. Unbalanced design is of course welcome.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. This is set to \eqn{10^{-9}} by default.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?  
}
\item{parallel}{
Do you want this to be executed in parallel or not. The parallel takes place in C++, and the number of threads 
is defined by each system's availiable cores.
}
%\item{ret}{
%This takes three values, 1, 2, or 3. The 1 returns the test statistic and the p-value. The 2 returns the BIC as 
%well and the 3 returns the BIC only. 
%}
\item{maxiters}{
The max number of iterations that can take place in each regression.
}
}

\details{
Many linear mixed models with a single covariate are fitted. We use Newton-Raphson as described in Demidenko (2013). The test statistic is the 
usual F-test. This model allows for random intercepts only.
}

\value{
A two-column matrix with the test statistics (Wald statistic) and the associated p-values 
(or their loggarithm). 
%If ret = 2 the matrix has three columns. The third one is the BIC.
%If ret = 3 the BIC is only returned.
}

\references{
Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons (excellent book). 
}

\author{
Stefanos Fafalios

R implementation and documentation: Stefanos Fafalios <stefanosfafalios@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{rint.reg}, \link{allbetas} \link{univglms}, \link{score.glms}, \link{logistic_only}
}
}

\examples{
## not a so good example
y <- rnorm(100)
id <- sample(1:10, 100, replace = TRUE)
x <- matrix( rnorm(100 * 100), ncol = 100)
a <- rint.regs(y, x, id)
x <- NULL
}
