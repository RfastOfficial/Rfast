\name{Limited number of eigenvalues and eigenvectors of a symmetric matrix}
\alias{eigen.sym}
\title{
Limited number of eigenvalues and eigenvectors of a symmetric matrix
}

\description{
Limited number of eigenvalues and eigenvectors of a symmetric matrix.
}

\usage{
eigen.sym(A, k, vectors = TRUE)
}

\arguments{
\item{A}{
A symmetric matrix.
}
\item{k}{
The number of eigenvalues and eigenvectors to extract.
}
\item{vectors}{
A flag that indicates if the eigenvectors will be returned (default: vectors = True)
}
}

\details{
The function calls the same function from the Armadillo library in C++. 
It is quite faster than R's built in function "eigen" if the number of 
eigenvalues and eigenvectors (argument k) is small. 

The k largest, in magnitude, eigenvalues are returned. Hence, if the matrix 
is not positive definite you may get negative eigenvalues as well. So, it is 
advised to use it with positive definite matrices.
}

\value{
A list including:
\item{values}{
The eigenvalues.
}
\item{vectors}{
The eigenvectors.
}
}

%\references{
%}

\author{
Armadillo library in C++ and Stefanos Fafalios and Manos Papadakis.

R implementation and documentation: Stefanos Fafalios <stefanosfafalios@gmail.com> and 
Manos Papadakis <papadakm95@gmail.com>. 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hd.eigen} 
}
}

\examples{
\donttest{
x <- matrnorm(500, 100 )
s <- Rfast::cova(x)
res<-eigen.sym(s, 5)
x <- s <- NULL
}
}


