\name{Bootstrap t-test for 2 independent samples}
\alias{boot.ttest2}
\title{
Bootstrap t-test for 2 independent samples
}

\description{
Bootstrap t-test for 2 independent samples.
}

\usage{
boot.ttest2(x, y, B = 999)
}

\arguments{
\item{x}{
A numerical vector with the data.
}
\item{y}{
A numerical vector with the data.
}
\item{B}{
The number of bootstrap samples to use.
}
}

\details{
Instead of sampling B times from each sample, we sample \eqn{sqrt{B}} from each of them and then take all pairs. 
Each bootstrap sample is independent of each other, hence there is no violation of the theory. 
}

\value{
A vector with the test statistic and the bootstrap p-value. 
}

\references{
B.L. Welch (1951). On the comparison of several mean values: an alternative approach. Biometrika, 38(3/4), 330-336.

Efron Bradley and Robert J. Tibshirani (1993). An introduction to the bootstrap. New York: Chapman \& Hall/CRC.

Chatzipantsiou C., Dimitriadis M., Papadakis M. and Tsagris M. (2019). 
Extremely efficient permutation and bootstrap hypothesis tests using R. 
To appear in the Journal of Modern Applied Statistical Methods.

https://arxiv.org/ftp/arxiv/papers/1806/1806.10947.pdf 
}

\author{
Michail Tsagris and Christina Chatzipantsiou

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Christina Chatzipantsiou <chatzipantsiou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ttest2}, \link{exact.ttest2}, \link{ftest}
}
}

\examples{
tic <- proc.time()
x <- rexp(40, 4)
y <- rbeta(50, 2.5, 7.5)
system.time( a <- boot.ttest2(x, y, 9999) )
a
}

