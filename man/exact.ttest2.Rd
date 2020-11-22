\name{Exact t-test for 2 independent samples}
\alias{exact.ttest2}
\title{
Exact t-test for 2 independent samples
}

\description{
Exact t-test for 2 independent samples.
}

\usage{
exact.ttest2(x, y)
}

\arguments{
\item{x}{
A numerical vector with the data.
}
\item{y}{
A numerical vector with the data.
}
}

\details{
This function performs an exact t-test. With few observations, permutation or bootstrap calculation
of the p-value is advisable. However, with even fewer observations, one can perform all possible 
permutations and calculate the exact p-value. This is what this function does. BUT, pay attention, 
as this works with few samples. If for example each sample contains 15 numbers, you will need a lot 
of memory (more than 17 GB) for this function to work. the reason is that we create the matrix with 
all possible permutations first and then perform the two-sample t-test. 
}

\value{
A vector with the number of permutations, test statistic and the permutation based p-value. 
}

\references{
B.L. Welch (1951). On the comparison of several mean values: an alternative approach. 
Biometrika, 38(3/4), 330-336.
}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and 
Manos Papadakis <papadakm95@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{boot.ttest2}, \link{ttest2}, \link{ftest}
}
}

\examples{
x <- rnorm(7)
y <- rnorm(7)
res<-exact.ttest2(x, y)
}

