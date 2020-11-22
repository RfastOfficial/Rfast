\name{Naive Bayes classifiers}
\alias{gaussian.nb}
\alias{poisson.nb}
\alias{multinom.nb}
\alias{geom.nb}
\alias{gammanb}

\title{
Naive Bayes classifiers
}

\description{
Gaussian, Poisson, geometric and multinomial naive Bayes classifiers.
}

\usage{
gaussian.nb(xnew = NULL, x, ina)
poisson.nb(xnew, x, ina)
multinom.nb(xnew, x, ina) 
geom.nb(xnew, x, ina, type = 1)
gammanb(xnew = NULL, x, ina, tol = 1e-07)   
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. For the Gaussian naive Bayes, 
this is set to NUUL, as you might want just the model and not to predict the membership of new observations. 
For the Gaussian case this contains any numbers, but for the multinomial and Poisson cases, the matrix must 
contain integer valued numbers only. 
}
\item{x}{
A numerical matrix with the observed predictor variable values. For the Gaussian case this contains any numbers, 
but for the multinomial and Poisson
cases, the matrix must contain integer valued numbers only. 
}
\item{ina}{
A numerical vector with strictly positive numbers, i.e. 1,2,3 indicating the groups of the dataset. 
Alternatively this can be a factor variable.
}
\item{type}{
Type 1 refers to the case where the minimum is zero and type 2 for the case of the minimum being 1. This 
is for the geometric distribution. This argument is for the geometric distribution. Type 1 refers to the case where the minimum is zero and type 2 for 
the case of the minimum being 1. 
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm in the gamma distribution.
}
}
%\details{
%
%}

\value{
For the Poisson and Multinomial naive Bayes classifiers the estimated group, a numerical vector with
1, 2, 3 and so on. For the Gaussian naive Bayes classifier a list including:
\item{mu}{
A matrix with the mean vector of each group based on the dataset.
}
\item{sigma}{
A matrix with the variance of each group and variable based on the dataset.
}
\item{ni}{
The sample size of each group in the dataset.
}
\item{est}{
The estimated group of the xnew observations. It returns a numerical value back regardless of the target 
variable being numerical as well or factor. Hence, it is suggested that you do \"as.numeric(target)\" in order to 
see what is the predicted class of the new data. 
}
For the Gamma classifier a list including:
\item{a}{
A matrix with the shape parameters.
}
\item{b}{
A matrix with the scale parameters.
}
\item{est}{
The estimated group of the xnew observations. It returns a numerical value back regardless of the target 
variable being numerical as well or factor. Hence, it is suggested that you do \"as.numeric(target)\" in order to 
see what is the predicted class of the new data. 
}
}

%\references{
%
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{gaussiannb.pred}, \link{colmeans}, \link{colVars} }
}

\examples{
x <- as.matrix(iris[, 1:4])
a <- gaussian.nb(x, x, iris[, 5])
x1 <- matrix( rpois(100 * 4, 5), ncol = 4)
x2 <- matrix( rpois(50 * 4, 10), ncol = 4)
x <- rbind(x1, x2)
ina <- c( rep(1, 100), rep(2, 50) )
res<-poisson.nb(x, x, ina)
res<-geom.nb(x, x, ina)
res<-multinom.nb(x, x, ina)
}
