\name{Many ANOVAS for count data with Poisson or quasi Poisson models}
\alias{colpoisson.anovas}
\alias{colquasipoisson.anovas}
\title{
Many ANOVAS for count data with Poisson or quasi Poisson models
}

\description{
Many ANOVAS for count data with Poisson or quasi Poisson models.
}

\usage{
colpoisson.anovas(y, x, logged = FALSE) 
colquasipoisson.anovas(y, x, logged = FALSE) 
}

\arguments{
\item{y}{
A numerical vector with the data.
}
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) 
and the columns are the variables. This must be a matrix with the categorical 
variables as numbers, starting from 1. Poisson or quassi Poisson ANOVA takes place
for each column.
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
}

\details{
Poisson or quassi Poisson ANOVA takes place at each column.
}

\value{
A matrix with the test statistic and the (logged) p-value for each predictor variable. 
In the case of the quasi Poisson, the \eqn{\phi} is returned as well. 
}

%\references{
%B.L. Welch (1951). On the comparison of several mean values: an alternative approach. 
%Biometrika, 38(3/4), 330-336.
%}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and 
Manos Papadakis <papadakm95@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{poisson.anova} \link{boot.ttest2}, \link{ttest2}, \link{ftest}
}
}

\examples{
y <- rpois(200, 10)
x <- matrix(rbinom(200 * 10, 3, 0.5 ), ncol = 10)
}

