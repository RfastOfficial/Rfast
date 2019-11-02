\name{Forward selection with generalised linear regression models}
\alias{fs.reg}
\title{
Variable selection in generalised linear regression models with forward selection
}

\description{
Variable selection in generalised linear regression models with forward selection
}

\usage{
fs.reg(y, ds, sig = 0.05, tol = 2, type = "logistic") 
}

\arguments{
\item{y}{
The dependent variable. This can either be a binary numeric (0, 1) or a vector with integers (numeric or integer class), count data. The first case is for
the binary logistic regression and the second for the Poisson regression. 
}
\item{ds}{
The dataset; provide a matrix where columns denote the variables and the rows the observations. The variables must be continuous, no categorical variables
are accepted. 
}
\item{sig}{
Significance level for assessing the p-values significance. Default value is 0.05.
}
\item{tol}{
The difference bewtween two successive values of the stopping rule. By default this is is set to 2. If for example, the BIC difference between two 
succesive models is less than 2, the process stops and the last variable, even though significant does not enter the model.
}
\item{type}{
If you have a binary dependent variable, put "logistic" or "quasibinomial". If you have percentages, values between 0 and 1, including 0 and or 1, use "quasibinomial" as well. 
If you have count data put "poisson".
}
}

\details{
The classical forward regression is implemented. The difference is that we have an extra step of check. Even if a variable is significant, the BIC of the model 
(with that variable) is calculated. If the decrease from the previous BIC (of the model without this variable) is less thatn a prespecified by the user value 
(default is 2) the variable wil enter. This way, we guard somehow against over-fitting. 
}

\value{
A matrix with for columns, the selected variables, the logarithm of their p-value, their test statistic and the BIC of the model with these variables included.
If no variable is selected, the matrix is empty.
}

%\references{
%Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. 
%In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
%}

\author{
Marios Dimitriadis

Documentation: Marios Dimitriadis <kmdimitriadis@gmail.com>.
}

\seealso{
\code{\link{cor.fsreg}, \link{logistic_only}, \link{poisson_only}, \link{glm_logistic}, \link{glm_poisson}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
x <- matrnorm(100, 50)
y <- rpois(100, 10)
a <- fs.reg(y, x, sig = 0.05, tol = 2, type = "poisson") 
x <- NULL
}

\keyword{ Forward regression }
\keyword{ variable selection }
\keyword{ binary data }
\keyword{ proportions }
\keyword{ Poisson }