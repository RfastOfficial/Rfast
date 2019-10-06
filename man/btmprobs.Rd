\name{Fitted probabilities of the Terry-Bradley model}
\alias{btmprobs}
\title{
Fitted probabilities of the Terry-Bradley model
}

\description{
Fitted probabilities of the Terry-Bradley model.
}

\usage{
btmprobs(x, tol = 1e-09)
}

\arguments{
\item{x}{
A numerical square, usually not symmetric, matrix with discrete valued data. Each entry is a frequency, to give an example, the number of wins. x[i, j] is the number 
of wins of home team i against guest team j. x[j, i] is the number of wins of home team j against guest team i. 
}
\item{tol}{
The tolerance level to terminate the iterative algorithm.
}
}

\details{
It fits a Bradley-Terry model to the given matrix and returns the fitted probabilities only.
}

\references{
Bradley R.A. and Terry M.E. (1952). Rank Analysis of Incomplete Block Designs: I. The Method of Paired Comparisons. Biometrika, 39(3/4):324-345.

Huang Tzu-Kuo, Ruby C. Weng and Chih-Jen Lin (2006). Generalized Bradley-Terry models and multi-class probability estimates. Journal of Machine Learning Research, 7:85-115.

Agresti A. (2002). Categorical Data Analysis (2nd ed). New York: Wiley.
}

\value{
A list including:
\item{iters}{
The numbetr of iterations required.
}
\item{probs}{
A vector with probabilities which sum to 1. This is the probability of win for each item (or team in our hypothetical example). 
}
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{g2tests}, \link{poisson.anova}, \link{anova}, \link{poisson_only}, \link{poisson.mle}
}
}

\examples{
x <- matrix( rpois(10 * 10, 10), ncol = 10) ## not the best example though
res<-btmprobs(x)
}

\keyword{ Bradley-Terry model }
\keyword{ fitted probabilities }
