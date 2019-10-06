\name{Prediction with some naive Bayes classifiers}
\alias{gaussiannb.pred}
\alias{poissonnb.pred}
\alias{multinomnb.pred}
\alias{gammanb.pred}
\alias{geomnb.pred}

\title{
Prediction with some naive Bayes classifiers
}

\description{
Prediction with some naive Bayes classifiers.
}

\usage{
gaussiannb.pred(xnew, m, s, ni)
poissonnb.pred(xnew, m)
multinomnb.pred(xnew, m)  
gammanb.pred(xnew, a, b)  
geomnb.pred(xnew, prob)  
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. For the Gaussian case this contains any numbers, but for the multinomial 
and Poisson cases, the matrix must contain integer valued numbers only. 
}
\item{m}{
A matrix with the group means. Each row corresponds to a group.
}
\item{s}{
A matrix with the group colum-wise variances. Each row corresponds to a group.
}
\item{ni}{
A vector with the frequencies of each group.
}
\item{a}{
A vector with the shape parameters of each group.
}
\item{b}{
A vector with the scale parameters of each group.
}
\item{prob}{
A vector with the sprobability parameters of each group.
}
}

%\details{
%
%}

\value{
A numerical vector with 1, 2, ... denoting the predicted group. 
}

%\references{
%
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{gaussian.nb}, \link{colpoisson.mle} \link{colVars} }
}

\examples{
ina <- sample(1:150, 100)
x <- as.matrix(iris[, 1:4])
id <- as.numeric(iris[, 5])
a <- gaussian.nb(xnew = NULL, x[ina, ], id[ina])
est <- gaussiannb.pred(x[-ina, ], a$mu, a$sigma, a$ni)
res<-table(id[-ina], est)
}

\keyword{ Discrimination }
\keyword{ naive Bayes }
\keyword{ normal distribution }
\keyword{ Poisson distribution }
\keyword{ multinomial distribution }