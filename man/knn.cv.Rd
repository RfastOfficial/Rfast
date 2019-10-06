\name{Cross-Validation for the k-NN algorithm}
\alias{knn.cv}
\title{
Cross-Validation for the k-NN algorithm
}

\description{
Cross-Validation for the k-NN algorithm.
}

\usage{
knn.cv(folds = NULL, nfolds = 10, stratified = FALSE, seed = FALSE, y, x, k, 
dist.type = "euclidean", type = "C", method = "average", freq.option = 0, 
pred.ret = FALSE, mem.eff = FALSE) 
}

\arguments{
\item{folds}{
A list with the indices of the folds.
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if "folds" is NULL.
}
\item{stratified}{
Do you want the folds to be selected using stratified random sampling? This preserves the analogy of the samples of each group. Make this 
TRUE if you wish, but only for the classification. 
If you have regression (type = "R"), do not put this to TRUE as it will cause problems or return wrong results.
}
\item{seed}{
If you set this to TRUE, the same folds will be created every time.
}
\item{y}{
A vector of data. The response variable, which can be either continuous or categorical (factor is acceptable).  
}
\item{x}{
A matrix with the available data, the predictor variables. 
}
\item{k}{
A vector with the possible numbers of nearest neighbours to be considered.
}
\item{dist.type}{
The type of distance to be used, "euclidean" or "manhattan".
}
\item{type}{
Do you want to do classification ("C") or regression ("R")? 
}
\item{method}{
If you do regression (type = "R"), then how should the predicted values be calculated? Choose among the average ("average"), median ("median") or 
the harmonic mean ("harmonic") of the closest neighbours.
}
\item{freq.option}{
If classification (type = "C") and ties occur in the prediction, more than one class have the same number 
of k nearest neighbours, there are three strategies available. Option 0 selects the first most frequent encountered. 
Option 1 randomly selects the most frequent value, in the case that there are duplicates. 
}
\item{pred.ret}{
If you want the predicted values returned set this to TRUE.
}
\item{mem.eff}{
Boolean value indicating a conservative or not use of memory. Lower usage of memory/Having this option on will lead to a slight 
decrease in execution speed and should ideally be on when the amount of memory in demand might be a concern.
}
}

\details{
The concept behind k-NN is simple. Suppose we have a matrix with predictor variables and a vector with the 
response variable (numerical or categorical). When a new vector with observations (predictor variables) is 
available, its corresponding response value, numerical or categorical, is to be predicted. Instead of using a 
model, parametric or not, one can use this ad hoc algorithm. 

The k smallest distances between the new predictor variables and the existing ones are calculated. In the 
case of regression, the average, median, or harmonic mean of the corresponding response values of these closest
predictor values are calculated. In the case of classification, i.e. categorical response value, a voting rule 
is applied. The most frequent group (response value) is where the new observation is to be allocated. 

This function does the cross-validation procedure to select the optimal k, the optimal number of nearest neighbours. 
The optimal in terms of some accuracy metric. For the classification it is the percentage of correct classification and for the
regression the mean squared error.
}

\value{
A list including:
\item{preds}{
If pred.ret is TRUE the predicted values for each fold are returned as elements in a list.
}
\item{crit}{
A vector whose length is equal to the number of k and is the accuracy metric for each k.
}
}

\references{
Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning. 
New York: Springer.

Cover TM and Hart PE (1967). Nearest neighbor pattern classification. IEEE Transactions on 
Information Theory. 13(1):21-27.

Tsagris Michail, Simon Preston and Andrew T.A. Wood (2016). Improved classification for compositional data using the 
\eqn{\alpha}-transformation. Journal of classification 33(2): 243-261.
}

\author{
Marios Dimitriadis

R implementation and documentation: Marios Dimitriadis <kmdimitriadis@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{knn}, \link{Dist}, \link{dista}, \link{dirknn.cv}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
y <- iris[, 5]
mod <- knn.cv(folds = NULL, nfolds = 10, stratified = TRUE, seed = FALSE, y = y, x = x, 
k = c(3, 4), dist.type = "euclidean", type = "C", method = "average", 
freq.option = 0, pred.ret = FALSE, mem.eff = FALSE) 
}

\keyword{ k-NN algorithm }
\keyword{ cross-validation }
