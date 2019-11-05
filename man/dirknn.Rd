\name{k-NN algorithm using the arc cosinus distance}
\alias{dirknn}
\title{
k-NN algorithm using the arc cosinus distance
}

\description{
It classifies new observations to some known groups via the k-NN algorithm.
}

\usage{
dirknn(xnew, x, y, k, type = "C", parallel = FALSE)
}

\arguments{
\item{xnew}{
The new data whose membership is to be predicted, a numeric matrix with unit vectors. 
In case you have one vector only make it a row vector (i.e. matrix with one row).
}
\item{x}{
The data, a numeric matrix with unit vectors.
}
\item{k}{
The number of nearest neighbours. It can also be a vector with many values.
}
\item{y}{
A numerical vector representing the class or label of each vector of x. 1, 2, 3, and so on.
It can also be a numerical vector with data in order to perform regression.
}
\item{type}{	
If your response variable y is numerical data, then this should be "R" (regression) or "WR"
for distance weighted based nearest neighbours. If y is in general categorical set this argument 
to "C" (classification) or to "WC" for distance weighted based nearest neighbours. 
}
\item{parallel}{
Do you want th ecalculations to take place in parallel? The default value is FALSE.
}
}

\details{
The standard algorithm is to keep the k nearest observations and see the groups of these observations. 
The new observation is allocated to the most frequent seen group. The non standard algorithm is to 
calculate the classical mean or the harmonic mean of the k nearest observations for each group. 
The new observation is allocated to the group with the smallest mean distance.

If you want regression, the predicted value is calculated as the average of the responses of the k nearest 
observations.
}
\value{
A matrix with the predicted group(s). It has as many columns as the values of k.
}

\author{
Stefanos Fafalios

R implementation and documentation: Stefanos Fafalios <stefanosfafalios@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dirknn.cv}, \link{knn}, \link{vmf.mle}, \link{spml.mle}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x/sqrt( rowSums(x^2) )
y<- as.numeric( iris[, 5] )
a <- dirknn(x, x, y, k = 2:10)
}

\keyword{ Directional k-NN algorithm }
\keyword{ Supervised classification }
