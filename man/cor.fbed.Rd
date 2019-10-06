\name{FBED variable selection method using the correlation}
\alias{cor.fbed}
\title{
FBED variable selection method using the correlation
}

\description{
FBED variable selection method using the correlation.
}

\usage{
cor.fbed(y, x, ystand = TRUE, xstand = TRUE, alpha = 0.05, K = 0)
}

\arguments{
\item{y}{
The response variable, a numeric vector.
}
\item{x}{
A matrix with the data, where the rows denote the samples and the columns are the variables. 
}
\item{ystand}{
If this is TRUE the response variable is centered. The mean is subtracted from every value.
}
\item{xstand}{
If this is TRUE the independent variables are standardised. 
}
\item{alpha}{
The significance level, set to 0.05 by default.
}
\item{K}{
The number of times to repeat the process. The default value is 0.
}
}

\details{
FBED stands for Forward Backward with Earcly Dropping. It is a variation of the classical forward selection, where at each step, only the 
statistically significant variables carry on. The rest are dropped. The process stops when no other variables can be selected. If K = 1, the process
is repeated testing sequentially again all those that have not been selected. If K > 1, then this is repeated. 

In the end, the backward selection is performed to remove any falsely included variables. This backward phase has not been implemented yet.  
}

\value{
A list including:
\item{runtime}{
The duration of the process.
}
\item{res}{
A matrix with the index of the selected variable, their test statistic value and the associated p-value.
}
\item{info}{
A matrix with two columns. The cumulative number of variables selected and the number of tests for each value of K.
}
}

\references{
Giorgos Borboudakis and Ioannis Tsamardinos (2017). Forward-Backward Selection with Early Dropping. Arxiv preprint: https://arxiv.org/pdf/1705.10770.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cor.fsreg}, \link{ompr}, \link{correls}, \link{fs.reg}
}
}

\examples{
x <- matrnorm(100, 100)
y <- rnorm(100)
a <- cor.fbed(y, x)
a
x <- NULL
}

