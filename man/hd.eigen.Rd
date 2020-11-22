\name{Eigenvalues and eigenvectors in high dimensional principal component analysis}
\alias{hd.eigen}
\title{
Eigenvalues in high dimensional principal component analysis
}

\description{
Eigenvalues in high dimensional (n<<p) principal component analysis.
}

\usage{
hd.eigen(x, center = TRUE, scale = FALSE, k = NULL, vectors = FALSE)
}

\arguments{
\item{x}{
A numerical \eqn{n \times p} matrix with data where the rows are the observations and the columns are the variables.
}
\item{center}{
Do you want your data centered? TRUE or FALSE.
}
\item{scale}{
Do you want each of your variables scaled, i.e. to have unit variance? TRUE or FALSE.
}
\item{k}{
If you want a specific number of eigenvalues and eigenvectors set it here, otherwise all 
eigenvalues (and eigenvectors if requested) will be returned.
}
\item{vectors}{
Do you want the eigenvectors be returned? By dafault this is FALSE.
}
}

\details{
When n<<p, at most the first n eigenvalues are non zero. Hence, there is no need to calculate the other 
p-n zero eigenvalues. When center is TRUE, the eigenvalues of the covariance matrix are calculated. 
When both the center and scale is TRUE the eigenvalues of the correlation matrix are calculated. One or more
eigenvectors (towards the end) will be 0. In general the signs might be the opposite than R's, but this makes
no difference. 
}

\value{
A list including:
\item{values}{
A vector with the n (or first k) eigenvalues. The divisor in the crossproduc matrix is n-1 and not n. 
}
\item{vectors}{
A matrix of \eqn{p \times n} or \eqn{p \times k} eigenvectors.
}
}

\author{
Michail Tsagris 

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{rmdp}
}
}

\examples{
x <- matrnorm( 40, 100)
a <- hd.eigen(x, FALSE, FALSE)
b <- prcomp(x, center = FALSE, scale = FALSE)
a
b$sdev^2
x <- NULL
}

\keyword{ Eigenvalues }
\keyword{ high dimensional data }
