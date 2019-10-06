\name{Binomial coefficient and its logarithm}
\alias{Lchoose}
\alias{Choose}
\title{
Binomial coefficient and its logarithm
}

\description{
Binomial coefficient and its logarithm.
}

\usage{
Lchoose(x, k)
Choose(x, k)
}

\arguments{
\item{x}{
A vector with integer values numbers.
}
\item{k}{
A positive non zero at most equal to x.
}
}

\details{
The binomial coefficient or its logarithm are evaluated.
}

\value{
A vector with the answers.
}


\author{
Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}


\seealso{
\code{ \link{comb_n}, \link{Lbeta}, \link{Lgamma} 
}
}

\examples{
x <- sample(20:30, 100, replace = TRUE)
res<-Choose(x, 4)
res<-Lchoose(x, 4)

x<-NULL
}

\keyword{ Factorials }
\keyword{ combinatorics}
