\name{Apply to each column a method under condition}
\alias{apply.condition}
\title{
Apply to each column a method under condition
}

\description{
Apply to each column a method under condition.
}

\usage{
apply.condition(x,method = "+",oper = ">",cond.val = 0)
}

\arguments{
\item{x}{
An integer matrix.
}
\item{method}{
One of: "+", "-", "*", "min", "max".
}
\item{oper}{
One of: ">, "<", ">=", "<=".
}
\item{cond.val}{
An integer value for the condition.
}
}

\details{
Apply to each col the specified method using the condition.
}

\value{
An integer vector with the coresponding values.
}


\author{
Manos Papadakis and Michail Tsagris

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> 
and Michail Tsagris <mtsagris@uoc.gr>.
}


\seealso{
\code{ \link{colsums}, \link{colMedians}, \link{colVars}
}
}

\examples{
x <- matrix(rpois(100,6),10, 10)
identical(apply(x,2,function(x){ sum(x[x>0]) }), apply.condition(x,"+",">",0))
x<-NULL
}


