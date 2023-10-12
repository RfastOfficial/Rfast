\name{Check Namespace and Rd files}
\alias{checkNamespace}
\alias{checkAliases}
\alias{checkExamples}
\alias{checkTF}
\alias{checkUsage}
\title{
Check Namespace and Rd files
}

\description{
Check Namespace/Rd and examples files.
}

\usage{
checkNamespace(path.namespace,path.rfolder)
checkAliases(path.man,path.rfolder)
checkTF(path.man)
checkExamples(path.man,package,each = 1,print.errors = stderr(),
	print.names = FALSE)
checkUsage(path.man,path.rfolder)
}

\arguments{
\item{path.namespace}{
An full path to the "NAMESPACE" file.
}
\item{package}{
A character vector with the name of the package.
}
\item{path.rfolder}{
An full path to the directory that contains the "R" files.
}
\item{path.man}{
An full path to the directory that contains the "Rd" files.
}
\item{each}{
An integer value for running \bold{each} example.
}
\item{print.errors}{
Print the errors to a file. By default it's "stderr()".
}
\item{print.names}{
A boolean value (TRUE/FALSE) for printing the names of the files before running the examples.
}
}

\details{
For function "checkNamespace": reads from the NAMESPACE folder all the export R functions, reads from folder R all the R functions and check if all the functions are export.

For function "checkAliases": reads from the man directory all the Rd files, then reads from each file the aliases and check if:
\itemize{
	\item All the R files has man file or an alias.
	\item All aliases belongs to functions.
	\item If there are dublicated aliases.
}

For function "checkExamples": reads from the man directory all the Rd files, then read from each file the examples and then run each of them. If you want to print the errors in any file then set "print.errors=file_name" or in the standard error "print.errors=stderr()" and then you will see all the errors for every file. For succeed run of your code you should first run "library(PACKAGE_NAME)". 
The argument "print.names" it is very helpful because if any of you function crashes R during running you will never know which one was. 
So setting it "TRUE", it will print the name of each file before runnign it's example.It might crash, but you will know which file. 
\bold{Remember that there is always an error timeout so it might didn't crash the current file but one from the previous.}

For function checkTF: reads from the man directory all the Rd files, then read from each file the examples and checks if any examples has the values "T" and "F" instead "TRUE" and "FALSE". The "T","F" is wrong.

For function checkUsage: reads from the man directory all the Rd files and for each man check if the usage section has the right signature for the functions from the R directory.

For functions "checkTF", "checkUsage", "checkAliases" you can choose which files not to read for both R and Rd. You must add in the first line of the file in comment the "attribute" "[dont read]". Then each function will know which file to read or not. For Rd you add "\%[dont read]" and for R "#[dont read]". Finally, these functions will return in the result a list of which files had this attribute.
}

\value{
For function "checkNamespace": a vector with the names of missing R files. (Don't use it for now)

For function "checkAliases": a list with 4 fields.
\itemize{
	\item \bold{Missing Man files}: A vector with the names of the missing Rd files or nothing.
	\item \bold{Missing R files}: A vector with the names of the missing R files or nothing.
	\item \bold{Duplicate alias}: A vector with the names of the dublicate aliases or nothing.
	\item{\bold{dont read}:}{ A list with 2 fields
		\itemize{
			\item R: A character vector whith the names of the files that had attribute "#[dont read]" or nothing.
			\item Rd: A character vector whith the names of the files that had attribute "\%[dont read]" or nothing.
		}
	}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

For function "checkExamples": a list with 3 fields 
\itemize{
	\item \bold{Errors}: A character vector with the names of the Rd files that produced an error.
	\item \bold{Big Examples}: A character vector with the names of the Rd files that has big examples per line.
	\item{\bold{dont read}:}{ A list with 2 fields
		\itemize{
			\item R: A character vector whith the names of the files that had attribute "#[dont read]" or nothing.
			\item Rd: A character vector whith the names of the files that had attribute "\%[dont read]" or nothing.
		}
	}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

For function "checkTF": a list with 3 fields 
\itemize{
	\item \bold{TRUE}: A character vector with the names of the Rd files that has "T" or nothing.
	\item \bold{FALSE}: A character vector with the names of the Rd files that has "F" or nothing.
	\item{\bold{dont read}:}{ A list with 2 fields
		\itemize{
			\item R: A character vector whith the names of the files that had attribute "#[dont read]" or nothing.
			\item Rd: A character vector whith the names of the files that had attribute "\%[dont read]" or nothing.
		}
	}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

For function "checkUsage": a list with 3 fields 
\itemize{
	\item \bold{missing functions}: A character vector with the name of the file that is missing and the Rd file that is found or nothing.
	\item \bold{missmatch functions}: A character vector with the name of the file that has missmatch function and the Rd file that is found or nothing.
	\item{\bold{dont read}:}{ A list with 2 fields
		\itemize{
			\item R: A character vector whith the names of the files that had attribute "#[dont read]" or nothing.
			\item Rd: A character vector whith the names of the files that had attribute "\%[dont read]" or nothing.
		}
	}
	\item \bold{hidden functions}: A character vector with the name of the functions thath have been declared as hidden.
	\item \bold{usage lines wider than 90 characters}: A list with the Rd's that have usages that are wider than 90 characters.
}
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{read.directory}, \link{AddToNamespace}, \link{sourceR}, \link{sourceRd}, \link{read.examples} 
}
}

\examples{
\donttest{
	#for example: path.namespace="C:\\some_file\\NAMESPACE"
	#for example: path.rfolder="C:\\some_file\\R\\"
	#for example: path.man="C:\\some_file\\man\\"
	#system.time( a<-checkNamespace(path.namespace,path.rfolder) )
	#system.time( b<-checkAliases(path.man,path.rfolder) )
	#system.time( b<-checkExamples(path.man) )
	#system.time( b<-checkExamples(path.man,2) )
	#system.time( b<-checkTF(path.man) )
	#system.time( b<-checkTF(path.man,path.rfolder) )
}
}

\keyword{ Namespace file }
\keyword{ Checking Alias }
\keyword{ Checking R }
\keyword{ Checking Rd }
\keyword{ Checking Examples }
\keyword{ Checking for TRUE }
\keyword{ Checking for FALSE }
\keyword{ Checking Usage section }