\name{Logistic and Poisson regression models}
\alias{glm_logistic}
\alias{glm_poisson}
\title{
Logistic and Poisson regression models
}

\description{
Logistic and Poisson regression models.
}

\usage{
glm_logistic(x, y, full = FALSE,tol = 1e-09, maxiters = 100)
glm_poisson(x, y, full = FALSE,tol = 1e-09)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns are the variables. 
This can be a matrix or a data.frame (with factors).
}
\item{y}{
The dependent variable; a numerical vector with two values (0 and 1) for the logistic regression
or integer values, 0, 1, 2,... for the Poisson regression. 
}
\item{full}{
If this is FALSE, the coefficients and the deviance will be returned only. If this is TRUE, more information is returned.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The max number of iterations that can take place in each regression.
}
}

\details{
The function is written in C++ and this is why it is very fast. 
}

\value{
When full is FALSE a list including:
\item{be}{
The regression coefficients.
}
\item{devi}{
The deviance of the model.
}
When full is TRUE a list including:
\item{info}{
The regression coefficients, their standard error, their Wald test statistic and their p-value.
}
\item{devi}{
The deviance.
}
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989.
}

\author{
Manos Papadakis <papadakm95@gmail.com>

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{poisson_only}, \link{logistic_only}, \link{univglms}, \link{regression}
}
}

\examples{
\dontrun{
x <- matrix(rnorm(100 * 3), ncol = 3)
y <- rbinom(100, 1, 0.6)   ## binary logistic regression
a1 <- glm_logistic(x, y, full = TRUE) 
a2 <- glm(y ~ x, binomial)

x <- matrix(rnorm(100 * 3), ncol = 3)
y <- rpois(100, 10)   ## binary logistic regression
b1 <- glm_poisson(x, y, full = TRUE) 
b2 <- glm(y ~ x, poisson)

x<-y<-a1<-a2<-b1<-b2<-NULL
}
}

\keyword{ Logistic regression }
\keyword{ poisson regression }
\keyword{ generalised linear models } 

