\name{Matrix multiplication}
\alias{mat.mult}
\alias{Crossprod}
\alias{Tcrossprod}

\title{
Matrix multiplication, Cross and Tcross product
}

\description{
Matrix multiplication, Cross and Tcross product.
}

\usage{
mat.mult(x, y)
Crossprod(x,y)
Tcrossprod(x,y)
}

\arguments{
\item{x}{
A numerical matrix.
}
\item{y}{
A numerical matrix.
}
}

\details{
The functions performs matrix multiplication, croos product and transpose cross product. There are faster(!) than R's function
for large matrices. Depending on the computer, maybe higher dimensions are 
required for the function to make a difference. The function runs in parallel
in C++. 
}

\value{
A matrix, the result of the matrix multiplication.
}

%\references{
%
%}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{transpose}, \link{colsums}
}
}

\examples{
\donttest{
x <- matrnorm(100, 100)
y <- matrnorm(100, 100)
a <- x %*% y
#b <- mat.mult(x, y) #Always execute in parallel
#b <- Crossprod(x, y) #Always execute in parallel
#b <- Tcrossprod(x, y) #Always execute in parallel
x <- NULL
y <- NULL
b <- NULL
}
}



