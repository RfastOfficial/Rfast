\name{Many univariate generalised linear models}
\alias{univglms}
\alias{univglms2}
\title{
Many univariate generalised linear regressions
}

\description{
It performs very many univariate generalised linear regressions.
}

\usage{
univglms(y, x, oiko = NULL, logged = FALSE)

univglms2(y, x, oiko = NULL, logged = FALSE)
}

\arguments{
\item{y}{
The dependent variable. It can be a factor or a numerical variable with two values only (binary logistic 
regression), a discrete valued vector (count data) corresponding to a poisson regression or a numerical 
vector with continuous values (normal regression). 
}
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the 
columns are the variables. For the "univglms" only continuous variables are allowed. You are advised to 
standardise the data before hand to avoid numerical overflow or similar issues. 
If you see NaN in the outcome, this might be the case. For the "univglms2" categorical variables are allowed
and hence this accepts data.frames. In this case, the categorical variables must be given as factor variables,
otherwise you might get wrong results.
}
\item{oiko}{
This can be either "normal", "poisson", "quasipoisson" or "binomial". If you are not sure leave it NULL and 
the function will check internally. However, you might have discrete data (e.g. years of age)
and want to perform many simple linear regressions. In this case you should specify the family. 
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
}

\details{
If you specify no family of distributions the function internally checkes the type of your data and 
decides on the type of regression to perform. The function is written in C++ and this is why it is very 
fast. It can accept thousands of predictor variables. It is usefull for univariate screening. We provide 
no p-value correction (such as fdr or q-values); this is up to the user. 
}

\value{
A matrix with the test statistic and the p-value for each predictor variable.
}

\references{
Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{logistic_only}, \link{poisson_only}, \link{allbetas}, \link{correls}, \link{regression}
}
}

\examples{
\dontrun{
x <- matrnorm(100, 50)
y <- rbinom(100, 1, 0.6)   ## binary logistic regression
a1 <- univglms(y, x) 
a2 <- glm(y ~ x[, 1], binomial)$deviance
a2 <- glm(y ~ 1, binomial)$null.dev - a2
x <- NULL
}
}



