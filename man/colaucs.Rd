\name{Many (and one) area under the curve values}
\alias{colaucs}
\alias{auc}
\title{
Many area under the curve values
}

\description{
Many area under the curve values.
}

\usage{
colaucs(group, preds)
auc(group, preds)
}

\arguments{
\item{group}{
A numerical vector with two values, one of which must be strictly 1. 
}
\item{preds}{
A numerical matrix with scores, probabilities or any other measure. 
In the case of auc this is a vector.
}
}

\details{
The AUCs are calculated column-wise or just an AUC if the vector function is used. 
}

\value{
A vector with length equal to the number of columns of the "preds" argument. 
The AUC vlaues for each column. If the "auc" function is used then a signle number is returned.
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and 
Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ttests}, \link{ttest}, \link{ftests}
}
}

\examples{
## 200 variables, hence 200 AUCs will be calculated
x <- matrix( rnorm(100 * 200), ncol = 200 )
ina <- rbinom(100, 1, 0.6)
system.time( colaucs(ina, x) )
a <- colaucs(ina, x) 
b <- auc(ina, x[, 1])
x <- NULL
}

\keyword{ Area aunder the curve}
