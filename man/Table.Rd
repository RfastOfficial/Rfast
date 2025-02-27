\name{Table Creation - Frequency of each value}
\alias{Table}
\alias{Table.sign}
\title{
Table Creation - Frequency of each value
}

\description{
Table Creation - Frequency of each value.
}

\usage{
Table(x,y=NULL,names = TRUE,useNA = FALSE,rm.zeros = FALSE)
Table.sign(x,names = TRUE,useNA = FALSE)
}

\arguments{
\item{x}{
A vector with numeric/character data.
}
\item{names}{
A logical value (TRUE/FALSE) for add names.
}
\item{y}{
A vector with numeric/character data. Doesn't work with "useNA".
}
\item{rm.zeros}{
A logical value for removing zero columns/rows. Only for integer vectors for now.
}
\item{useNA}{
Table: 
	Integer/logical value:

	 FALSE: not NA values in vector.
	 TRUE: count NAs and add the value in the last position of the returned vector.
	 any other integer except 0,1: for just removing NAs.

Table.sign:
	Logical value, TRUE, for count NAs. Otherwise FALSE.

Doesn't work character data.
}
}

\details{
Like R's "table":

for giving one argument,"x":
	 If "names" is FALSE then, if "useNA" is TRUE  then the NAs will be count, if is  FALSE it means there are no NAs and for any other integer value the NAs will be ignored.

for giving two arguments,"x","y":
	If "names" is FALSE then, creates the contigency table, otherwise sets the col-row names with discrete values. If "rm.zeros" is FALSE then it won't remove the zero columns/rows from the result but it will work only for positive integers for now. For this if "names" is TRUE then the col-row names will be the seq(min(),max()) for "x","y". In future updates it will be changed. 

for both algorithms: 
	You can't use "useNA" with "names" for now. 
	It is much faster to get the result without names (names = FALSE) but all the algorithms are more efficient than R's. 

Like R's "table(sign())" but more efficient. Count the frequencies of positives, negatives, zeros and NAs values. If argument "names" is FALSE then the returned vector doesn't have names. Otherwise "-1,0,+1,NA".
If "useNA" is TRUE then the NAs will be count, otherwise not. You can use "useNA" with "names".
}

\value{
Table: 

for giving one argument,"x":	
	if "names" is TRUE then return a vector with names the discrete values of "x" and values there frequencies, otherwise only the frequencies

for giving two arguments,"x","y":
	if "names" is TRUE then return a contigency matrix with rownames the discrete values of "x", colnames the dicrete values of "y" and values the freuquencies of the pairs, otherwise only the freuquencies of the pairs.

Table.sign:
	A vector with 4 values/frequencies:
		index 1: negatives
		index 2: zeros
		index 3: postives
		if "names" is TRUE then the returned vector have names "-1,0,+1".
		if "useNA" is TRUE then 4th value has the frequencies of NAs and the returned vector will have one more name, "-1,0,+1,NA", if "names" is also TRUE.


}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}


\seealso{
\code{\link{colShuffle}, \link{colVars}, \link{colmeans}, \link{is_integer}, \link{as_integer}
}
}

\examples{
x<-runif(10)
y1<-Table(x)
y2<-as.vector(table(x)) # Neads a lot of time.
all.equal(y1,y2)
y1<-Table(x,names=FALSE)
all.equal(y1,y2) # the name attribute of y1 is null 

y1<-Table.sign(x)
y2<-table(sign(x))
all.equal(y1,y2)

x<-y1<-y2<-NULL
}

\keyword{ Table Creation }