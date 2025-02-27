\name{Source many R files}
\alias{sourceR}
\alias{sourceRd}
\title{
Source many R files
}

\description{
Source many R/Rd files.
}

\usage{
sourceR(path,local=FALSE,encode = "UTF-8",print.errors=FALSE)
sourceRd(path,print.errors=FALSE, macros = NULL)
}

\arguments{
\item{path}{
An full path to the directory where R file are.
}
\item{local}{
TRUE, FALSE or an environment, determining where the parsed expressions are evaluated. FALSE (the default) corresponds to the user's workspace (the global environment) and TRUE to the environment from which source is called.
}
\item{encode}{
Character vector. The encoding(s) to be assumed when file is a character string: see file. A possible value is "unknown" when the encoding is guessed: see the "Encodings" section.
}
\item{print.errors}{
A boolean value (TRUE/FALSE) for printing the errors, if exists, for every file.
}
\item{macros}{
A full or relative path to the macros file.
}
}

\details{
Reads many R files and source them.
}

\value{
Returns the files that had produced errors during source.
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\link{AddToNamespace} 
}

\examples{
# for example: path="C:\\some_file\\R\" where is R files are
# system.time( a<-sourceR(path) )
# for example: path="C:\\some_file\\man\" where is Rd files are
# system.time( a<-sourceRd(path) )
}

\keyword{ Namespace file }
\keyword{ Export functions }
