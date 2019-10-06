\name{Backward selection regression}
\alias{bs.reg}
\title{
Backward selection regression
}

\description{
Backward selection regression.
}

\usage{
bs.reg(y, x, alpha = 0.05, type = "logistic")
}

\arguments{
\item{y}{
A numerical vector with the response variable values. It can either be of 0 and 1 values (Logistic regression) or of integer values 0, 1, 2,... (Poisson regression).
}
\item{x}{
A numerical matrix with the candidate variables.
}
\item{alpha}{
Threshold (suitable values are in [0,1]) for assessing the significance of p-values. The default value is at 0.05.
}
\item{type}{
For the Logistic regression put "logistic" (default value) and for Poisson type "poisson".
}
}

\details{
This function currently implements only the binary Logistic and Poisson regressions. If the sample size is less than the number of variables a notification message 
will appear and no backward regression will be performed. 
}

\value{
The output of the algorithm is an S3 object including:
\item{info}{
A matrix with the non selected variables and their latest test statistics and p-values.
}
\item{Vars}{
A vector with the selected variables.
}
}

\author{
Marios Dimitriadis
 
R implementation and documentation: Marios Dimitriadis <mtsagris@csd.uoc.gr>
}

\seealso{
\code{\link{fs.reg}, \link{univglms}, \link{cor.fsreg} }
}

\examples{
y <- rbinom(50, 1, 0.5)
x <- matrnorm(50, 10)
res<-bs.reg(y, x)
}

