\name{Many two-way ANOVAs}
\alias{twoway.anovas}
\title{
Many two-way ANOVAs
}

\description{
Many two-way ANOVAs.
}

\usage{
twoway.anovas(y, x1, x2, interact = FALSE, logged = FALSE)
}

\arguments{
\item{y}{
A matrix with the data, where the rows denote the observations (and the two groups) and the columns are the variables.
}
\item{x1}{
A numerical vector with 1s, 2s, 3s and so one indicating the two groups. Alternatively it can be a factor variable. 
This is the one factor.
}
\item{x2}{
A numerical vector with 1s, 2s, 3s and so one indicating the two groups. Alternatively it can be a factor variable. 
This is the other factor.
}
\item{interact}{
A boolean variable specifying whether you want to test for interaction.  
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
The classical two-way ANOVA design is performed. Note that the design must be balanced. For every combination
of values of the two factors, x1 and x2 the same number of observations must exist. If that's not the case, 
regression models must be used. 
}

\value{
A matrix with the test statistic and the p-value of each test. 
}

\references{
D.C. Montgomery (2001). Design and analysis of experiments (5th Edition). New York: John Wiley \& Sons 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ancovas}, \link{ftests}, \link{ttests}
}
}

\examples{
y <- as.matrix( rnorm(125) )
x1 <- rep(1:5, 25)
x2 <- rep(1:5, each = 25)
x1 <- factor(x1)
x2 <- factor(x2)  
res<-anova( lm(y ~ x1 + x2) )
res<-twoway.anovas(y, x1, x2)
res<-anova( lm(y ~ x1*x2) )
res<-twoway.anovas(y, x1, x2, interact = TRUE) 
y <- matrnorm(125, 100)
system.time( a1 <- twoway.anovas(y, x1, x2) )
system.time( a2 <- twoway.anovas(y, x1, x2, interact = TRUE) )
y <- NULL
}

\keyword{ Two-way ANOVA }
\keyword{ interaction }
