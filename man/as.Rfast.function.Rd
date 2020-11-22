\name{Convert R function to the Rfast's coresponding}
\alias{as.Rfast.function}
\title{
Convert R function to the Rfast's coresponding
}

\description{
Convert R function to the Rfast's coresponding.
}

\usage{
as.Rfast.function(Rfunction.name,margin=NULL)
}

\arguments{
\item{Rfunction.name}{
An character value with the name of the function.
}
\item{margin}{
A logical function for return the column-row wise function.
}
}

\details{
Given the name of R function, it returns the coresponding function's name from Rfast.
}

\value{
The coresponding Rfast function.
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

res<-as.Rfast.function("var")

}


