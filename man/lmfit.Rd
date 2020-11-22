\name{Linear models for large scale data}
\alias{lmfit}
\title{
Linear models for large scale data
}

\description{
Linear models for large scale data.
}

\usage{
lmfit(x, y, w = NULL)
}

\arguments{
\item{x}{
The design matrix with the data, where each column refers to a different sample of subjects. 
You must supply the design matrix, with the column of 1s. This function is the analogue of 
lm.fit and .lm.fit. 
}
\item{y}{
A numerical vector or a numerical matrix. 
}
\item{w}{
An optional numerical vector with weights. Note that if you supply this, the function does not make them sum to 1. So, you should do it. 
}
}

\details{
We have simply exploitted R's powerful function and managed to do better than .lm.fit which is a really powerful 
function as well. This is a bare bones function as it returns only two things, the coefficients and the residuals. 
\link{.lm.fit} returns more and \link{lm.fit} even more and finally \link{lm} returns too much.  
The motivatrion came form this site https://m-clark.github.io/docs/fastr.html . We changed the function a bit.
}

\value{
A list including:
\item{be}{
The beta coefficients.
}
\item{residuals}{
The residuals of the linear model(s). 
}
}

\references{
Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{regression}, \link{allbetas}, \link{correls}, \link{mvbetas}, \link{cor.fsreg}
}
}

\examples{
n <- 200  ;  p <- 5
X <- matrnorm(n, p)
y <- rnorm(200)
a1 <- .lm.fit(X, y) 
a2 <- lmfit(X, y) 
x <- NULL
}

\keyword{ Linear models }
\keyword{ large scale data}  
