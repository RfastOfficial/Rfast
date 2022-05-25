\name{ Insert/remove function names in/from the NAMESPACE file }
\alias{ AddToNamespace }
\alias{ RemoveFromNamespace }
\title{
Insert/remove function names in/from the NAMESPACE file
}

\description{
Insert/remove function names in/from the NAMESPACE file.
}

\usage{
AddToNamespace(path.namespace,path.rfolder)
RemoveFromNamespace(path.namespace,files.to.remove)
}

\arguments{
\item{path.namespace}{
An full path to the NAMESPACE file.
}
\item{path.rfolder}{
An full path to the directory the new files to be added are stored.
}
\item{files.to.remove}{
An character with the names of the functions to be removed from file NAMESPACE.
}
}

\details{
AddToNameSpace: Reads the files that are exported in NAMESPACE and the functions that are inside rfolder (where R files are) and insert every function that is not exported. For that you must add the attribute "#[export]" above every function you wish to export. Also you can use the attribute "#[export s3]" for exporting S3methods. Finally, if you don't want the program to read a file just add at the top of the file the attribute "#[dont read]".

RemoveFromNamespace: Remove every function, from argument "files.to.remove", from NAMESPACE.
}

\value{
AddToNameSpace:
\item{without export}{
A character vector with the names of the R functions that don't have te "#[export]" attribute.
}
\item{hidden functions}{
A character vector with the names of the R functions that are hidden.
}

RemoveFromNamespace: Return the files that could not be removed.
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colShuffle}, \link{colVars}, \link{colmeans}, \link{read.directory}
}
}

\examples{
\dontrun{
	for example: path.namespace="C:\\some_file\\NAMESPACE" where is NAMESPACE file
	path.rfolder="C:\\some_file\\R\\" where is R files are
	system.time( a<-AddToNamespace(path.namespace,path.rfolder) )
	if(length(a)==0){
		print("all the files are inserted")
	}else{
		print("The new files that inserted are: \n")
		a
	}
	system.time( a<-RemoveFromNamespace(path.namespace,c("a","b")) )
	if(length(a)==0){
		print("all the files are inserted")
	}else{
		print("The files thatcould not be deleted are: \n")
		a
	}
}
}

\keyword{ Namespace file }
\keyword{ Export functions }
\keyword{ Remove functions }
