\name{Empirical and exponential empirical likelihood tests for one sample}
\alias{eel.test1}
\alias{el.test1}
\title{
Empirical and exponential empirical likelihood tests for one sample
}

\description{
Empirical and exponential empirical likelihood tests for one sample.
}

\usage{
eel.test1(x, mu, tol = 1e-09, logged = FALSE)
el.test1(x, mu, tol = 1e-07, logged = FALSE) 
}

\arguments{
\item{x}{
A numerical vector. 
}
\item{mu}{
The hypothesised mean value.
}
\item{tol}{
The tolerance value to stop the iterations of the Newton-Raphson.
}
\item{logged}{
Should the logarithm of the p-value be returned? TRUE or FALSE.
}
}

\details{
Exponential empirical likelihood is a non parametric method. In this case we use it as the non parametric alternative to the t-test.
Newton-Raphson is used to maximise the log-likelihood ratio test statistic. In the case of no solution, NULL is returned. 
Despite the function having beeen written in R, it is pretty fast. As for the empirical likelihood ratio test, there is a condition for the 
range of possible values of mu. If mu is outside this range it is rejected immediately. 
}

\value{
\item{iters}{
The number of iterations required by the Newton-Raphson algorithm. If no covnergence occured this is NULL. This is not returned for the
empircial likelihood ratio test.
}
\item{info}{
A vector with three elements, the value of the \eqn{\lambda}, the likelihood ratio test statistic and the relevant p-value.
If no convergence occured, the value of the \eqn{\lambda} before is becomes NA, the value of test statistic is \eqn{10^5} and the p-value is 0.
No convergence can be interpreted as rejection of the hypothesis test.
}
\item{p}{
The estimated probabilities, one for each observation. If no covnergence occured this is NULL.
}
}

\references{
Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ftest}, \link{ttest1}
}
}

\examples{
x <- rnorm(500)
system.time(a1 <- eel.test1(x, 0) )
system.time(a2 <- el.test1(x, 0) )
}

\keyword{ Hypothesis testing }
\keyword{ one sample }
\keyword{ non parametric test }
