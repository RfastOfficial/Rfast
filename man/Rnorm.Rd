\name{Simulation of random values from a normal distribution}
\alias{Rnorm}
\title{
Simulation of random values from a normal distribution
}

\description{
Simulation of random values from a normal distribution.
}

\usage{
Rnorm(n, m = 0, s = 1, seed = NULL)
}

\arguments{
\item{n}{
The sample size. 
}
\item{m}{
The mean, set to 0 by default.
}
\item{s}{
The standard devation, set to 1 by default.
}
\item{seed}{
If you want the same to be generated again use a seed for the generator, an integer number.
}
}

\details{
By using the Ziggurat method of generating standard normal variates, this function is really fast when you 
want to generate large vectors. For less than 2,000 this might make no difference when compared with R's 
"rnorm", but for 10,000 this will be 6-7 times faster.
}

\value{
A vector with n values. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{matrnorm}, \link{rvonmises}, \link{rvmf}, \link{rmvnorm}
}
}

\examples{
x <- Rnorm(500)
}




