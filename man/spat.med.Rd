\name{Spatial median for Euclidean data}
\alias{spat.med}
\title{
Spatial median for Euclidean data
}

\description{
Spatial median for Euclidean data.
}

\usage{
spat.med(x, tol = 1e-09)
}

\arguments{
\item{x}{
A matrix with Euclidean data, continuous variables.
}
\item{tol}{
A tolerance level to terminate the process. This is set to 1e-09 by default.
}
}

\details{
The spatial median, using a fixed point iterative algorithm, for Euclidean data is calculated. It is a robust location estimate.
}

\value{
A vector with the spatial median.
}

\references{
Jyrki Mottonen, Klaus Nordhausen and Hannu Oja (2010). Asymptotic theory of the spatial median.
In Nonparametrics and Robustness in Modern Statistical Inference and Time Series Analysis:
A Festschrift in honor of Professor Jana Jureckova.

T. Karkkaminen and S. Ayramo (2005). On computation of spatial median for robust data mining.
Evolutionary and Deterministic Methods for Design, Optimization and Control with Applications to
Industrial and Societal Problems, EUROGEN 2005, R. Schilling, W.Haase, J. Periaux, H. Baier, G. Bugeda (Eds)
FLM, Munich. http://users.jyu.fi/~samiayr/pdf/ayramo_eurogen05.pdf
}

\author{
Manos Papadakis and Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{colMedians}
}
}

\examples{
res<-spat.med( as.matrix( iris[, 1:4] ) )
res<-colMeans( as.matrix(iris[, 1:4]) )
res<-colMedians( as.matrix(iris[, 1:4]) )
}

\keyword{ spatial median }
\keyword{ robust statistics }
