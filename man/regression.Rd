\name{Many univariate simple linear regressions}
\alias{regression}
\title{
Many univariate simple linear regressions
}

\description{
It performs very many univariate simple linear regressions with or without categorical variables.
}

\usage{
regression(x, y, poia = NULL, logged = FALSE)
}

\arguments{
\item{x}{
A data.frame or a matrix with the data, where the rows denote the samples (and the two groups) and the columns 
are the variables. A data frame is expected if you have categorical predictor variables. If you only have continuous 
predictor variables you should the function \code{\link{allbetas}} instead as it is faster. 
}
\item{y}{
The dependent variable; a numerical vector.
}
\item{poia}{
If the "x" is a data.frame and you know the indices of the columns which are categorical variables supply it here.
}
\item{logged}{
Do you want the logarithm of the p-values be returned? The default value is FALSE.
}
}

\details{
Some parts of the function will be transferred in C++. It can accept thousands of predictor variables. 
It is usefull for univariate screening. We provide no p-value correction (such as fdr or q-values); this is 
up to the user. 
}

\value{
A matrix with two columns, the test statistic value and its corresponding (logged) p-value.
}

\references{
Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Manos Papadakis <papadakm95@gmail.com>

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{univglms}, \link{allbetas}, \link{correls}, \link{univglms}, \link{mvbetas}
}
}

\examples{
y <- rnorm(150)   
a <- regression(iris, y)
a
summary(lm(y ~ iris[, 5]) )  ## check the F-test
}

\keyword{ Simple linear regressions }
\keyword{ categorical variables }

