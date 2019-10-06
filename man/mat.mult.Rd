\name{Matrix multiplication}
\alias{mat.mult}

\title{
Matrix multiplication
}

\description{
Matrix multiplication.
}

\usage{
mat.mult(x, y)

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
The function performs matrix multiplication. It is faster(!) than R's function
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
Manos Papadakis and Michail Tsagris

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com> 
and Michail Tsagris <mtsagris@yahoo.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{transpose}, \link{colsums}
}
}

\examples{
\dontrun{
x <- matrnorm(100, 100)
y <- matrnorm(100, 100)
a <- x %*% y
b <- mat.mult(x, y)
x <- NULL
y <- NULL
}
}



