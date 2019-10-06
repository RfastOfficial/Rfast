\name{k nearest neighbours algorithm (k-NN)}
\alias{knn}
\title{
k nearest neighbours algorithm (k-NN)
}

\description{
k nearest neighbours algorithm (k-NN). 
}

\usage{
knn(xnew, y, x, k, dist.type = "euclidean", type = "C", method = "average", 
freq.option = 0, mem.eff = FALSE)
}

\arguments{
\item{xnew}{
The new data, new predictor variable values. A matrix with numerical data. 
}
\item{y}{
A vector with the response variable, whose values for the new data we wish to predict. This can be 
numerical data, factor or discrete, 0, 1, ... The latter two cases are for classification.
}
\item{x}{
The dataset. A matrix with numerical data.
}
\item{k}{
The number of nearest neighbours to use. The number can either be a single value or a vector with multiple values.
}
\item{dist.type}{
The type of distance to be used. Either \"euclidean\" or \"manhattan\". 
}
\item{type}{
If your response variable \"y\" is numerical data, then this should be \"R\" (regression). If \"y\" is in general 
categorical, factor or discrete set this argument to \"C\" (classification).
}
\item{method}{
In case you have regression (type = \"R\") you want a way to summarise the prediction. If you want to take the 
average of the reponses of the k closest observations, type \"average\". For the median, type \"median\" and for the 
harmonic mean, type \"harmonic\".
}
\item{freq.option}{
If classification (type = \"C\") and ties occur in the prediction, more than one class has the same number 
of k nearest neighbours, in which case there are two strategies available: Option 0 selects the first most frequent encountered. Option 1 randomly selects the most frequent value, in the case that there are duplicates. 
}
\item{mem.eff}{
Boolean value indicating a conservative or not use of memory. Lower usage of memory/Having this option on will lead to a slight 
decrease in execution speed and should ideally be on when the amount of memory in demand might be a concern.
}
}

\details{
The concept behind k-NN is simple. Suppose we have a matrix with predictor variables and a vector with the 
response variable (numerical or categorical). When a new vector with observations (predictor variables) is 
available, its corresponding response value, numerical or category is to be predicted. Instead of using a 
model, parametric or not, one can use this ad hoc algorithm. 

The k smallest distances between the new predictor variables and the existing ones are calculated. In the 
case of regression, the average, median or harmonic mean of the corresponding respone values of these closest
predictor values are calculated. In the case of classification, i.e. categorical response value, a voting rule 
is applied. The most frequent group (response value) is where the new observation is to be allocated. 
}

\value{
A matrix whose number of columns is equal to the size of k. If in the input you provided there is just one value of
k, then a matrix with one column is returned containing the predicted values. If more than one value was supplied,
the matrix will contain the predicted values for every value of k. 
}

\references{
Cover TM and Hart PE (1967). Nearest neighbor pattern classification. IEEE Transactions on 
Information Theory. 13(1):21-27.

Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning. 
New York: Springer.

http://web.stanford.edu/~hastie/ElemStatLearn/printings/ESLII_print12.pdf 	

http://statlink.tripod.com/id3.html
}

\author{
Marios Dimitriadis

R implementation and documentation: Marios Dimitriadis <kmdimitriadis@gmail.com>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{knn.cv}, \link{dirknn}, \link{logistic_only}, \link{fs.reg}, \link{cor.fsreg} }
}

\examples{
# Simulate a dataset with continuous data
x <- as.matrix(iris[, 1:4])
y <- as.numeric(iris[, 5])
id <- sample(1:150, 120)
mod <- knn(x[-id, ], y[id], x[id, ], k = c(4, 5, 6), type = "C", mem.eff = FALSE)
mod # Predicted values of y for 3 values of k. 
res<-table(mod[, 1], y[-id] )  # Confusion matrix for k = 4
res<-table(mod[, 2], y[-id] )  # Confusion matrix for k = 5
res<-table(mod[, 3], y[-id] )  # Confusion matrix for k = 6
}

\keyword{ k-NN algorithm }
