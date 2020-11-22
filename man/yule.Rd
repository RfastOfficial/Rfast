\name{Yule's Y (coefficient of colligation)}
\alias{yule}
\title{
Yule's Y (coefficient of colligation)
}

\description{
Yule's Y (coefficient of colligation).
}

\usage{
yule(x)
}

\arguments{
\item{x}{
A 2 x 2 matrix or a vector with 4 elements. In the case of the vector make sure it corresponds to the correct table.
}
}

\details{
Yule's coefficient of colligation is calculated.
}

\value{
Yule's Y is returned.
}

\references{
Yule G. Udny (1912). On the Methods of Measuring Association Between Two Attributes. 
Journal of the Royal Statistical Society, 75(6):579-652.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{col.yule}, \link{odds.ratio}
}
}

\examples{
x <- rpois(4, 30) + 2
res<-yule(x)
res<-yule( matrix(x, ncol = 2) )
}


