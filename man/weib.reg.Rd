\name{Weibull regression model}
\alias{weib.reg}
\title{
Weibull regression model
}

\description{
Weibull regression model.
}

\usage{
weib.reg(y, x, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable; a numerical vector with strictly positive data, i.e. greater than zero. 
}
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns are the variables. 
This can be a matrix or a data.frame (with factors).
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The max number of iterations that can take place in each regression.
}
}

\details{
The function is written in C++ and this is why it is very fast. No standard errors are returned as they are not 
corectly estimated. We focused on speed.
}

\value{
When full is FALSE a list including:
\item{iters}{
The iterations required by the Newton-Raphson.
}
\item{loglik}{
The log-likelihood of the model.
}
\item{shape}{
The shape parameter of the Weibull regression.
}
\item{be}{
The regression coefficients. 
}
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989.
}

\author{
Stefanos Fafalios

R implementation and documentation: Stefanos Fafalios <stefanosfafalios@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{poisson_only}, \link{logistic_only}, \link{univglms}, \link{regression}
}
}

\examples{
\donttest{
x <- matrix(rnorm(100 * 2), ncol = 2)
y <- rexp(100, 1)
res<-weib.reg(y, x)
}
}

