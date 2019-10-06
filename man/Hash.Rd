\name{Hash object}
\alias{Hash}
\alias{Hash.key.multi}
\alias{print.Hash}
\alias{length.Hash}
\alias{[.Hash}
\alias{[<-.Hash}
\title{
Hash object
}

\description{
Hash object.
}

\usage{
Hash(keys=NULL,values=NULL)
Hash.key.multi(x,...,sep = " ")
\method{[}{Hash}(x,...,sep = " ") <- value
\method{[}{Hash}(x,...,sep = " ")
\method{print}{Hash}(x,...)
\method{length}{Hash}(x)
}

\arguments{
\item{x}{
A Hash object, using Hash function.
}
\item{values}{
A vector with the values you want to store.
}
\item{value}{
The values you want to store.
}
\item{keys}{
A vector with keys for each values.
}
\item{sep}{
A character value using to separate the multiple keys for each value.
}
\item{...}{
One or more values for access or find elements.
}
}

\details{
If you want to delete a key just insert the global variable "Rfast:::delete".

Hash: Create Hash object where every key has a value. Specify the type from the beggining (for speed). Use the argument "type" with one of the values "new.env, logical, character, integer, numeric".
Hash.key.multi: search if key exists. If the keys are multiple, then use the argument "substr" to search inside each multiple for the specific key.
}

\value{
A Hash object.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{  \link{hash.list}, \link{hash.find} }
}

\examples{
x <- Hash(rnorm(10),sample(1:10))

x[1,2,13] <- 0.1234 # insert value using multi key. the same as x["1 2 13"] <- 0.1234
x[1,2,3] <- 15 # insert value using multi key. the same as x["1 2 3"] <- 15

Hash.key.multi(x,"1")
x # print Hash object using S3 generic
#x[1,2,3] <- Rfast:::delete # delete multi key. the same as x["1 2 3"] <- NULL
length(x)
}

\keyword{ Hash tables }
\keyword{ list }
