\name{Fast and general - untyped represantation of a factor variable}
\alias{ufactor}
\alias{print.ufactor}
\alias{[.ufactor}
\title{
Fast and general represantation of a factor variable
}

\description{
Fast and general represantation of a factor variable.
}

\usage{
ufactor(x)
\method{[}{ufactor}(x,i)
\method{print}{ufactor}(x,...)
}

\arguments{
\item{x}{
A vector with data. 
}
\item{i}{
An integer value/vector which is the index/indices to the element you want to access.
}
\item{...}{
Anything the user want.
}
}

\details{
This is a general implementation of factor structure. For access the fields of a "ufactor" use the "$" operator.
}

\value{
An object of class "ufactor". This object holds 2 fields:

levels: the levels of the variable in his initial type
values: the values of the variable in his initial type
}


\author{
Manos Papadakis

R implementation and documentation: and Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colVars}, \link{factor}
}
}

\examples{
x <- rnorm(10)
R.factor<- as.factor(x)
Rfast.factor <- ufactor(x)

identical(levels(R.factor),Rfast.factor$levels) # TRUE
identical(as.numeric(R.factor),Rfast.factor$values) # TRUE
x<-R.factor<-Rfast.factor<-NULL
}

\keyword{ Factor }
