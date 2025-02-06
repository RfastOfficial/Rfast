\name{Permutation based p-value for the Pearson correlation coefficient}
\alias{permcor}
\title{
Permutation based p-value for the Pearson correlation coefficient
}

\description{
Permutation based p-value for the Pearson correlation coefficient.
}

\usage{
permcor(x, y, R = 999) 
}

\arguments{
\item{x}{
A numerical vector with the first variable. 
}
\item{y}{
A numerical vector with the second variable. 
}
\item{R}{
The number of permutations to be conducted; set to 999 by default. 
}
}

\details{
This is a very low computational calculation of the p-value. Try it yourselves. 
}

\value{
A vector consisting of two values, the Pearson correlation and the permutation based p-value. 
}

\references{
Chatzipantsiou C., Dimitriadis M., Papadakis M. and Tsagris M. (2019). 
Extremely efficient permutation and bootstrap hypothesis tests using R. 
To appear in the Journal of Modern Applied Statistical Methods.

https://arxiv.org/ftp/arxiv/papers/1806/1806.10947.pdf
}

\author{
Marios Dimitriadis and Michail Tsagris

R implementation and documentation: Marios Dimitriadis and Michail Tsagris <kmdimitriadis@gmail.com> 
and <mtsagris@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pc.skel}}
}

\examples{
x <- iris[, 1]
y <- iris[, 2]
permcor(x, y, R = 499)
}
