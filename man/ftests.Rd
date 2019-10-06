\name{Many multi-sample tests}
\alias{ftests}
\alias{anovas}
\alias{vartests}
\alias{block.anovas}
\title{
Many multi-sample tests
}

\description{
Many multi-sample tests.
}

\usage{
ftests(x, ina, logged = FALSE)
anovas(x, ina, logged = FALSE)
vartests(x, ina, type = "levene", logged = FALSE)
block.anovas(x, treat, block, logged = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the observations (and the two groups) and the columns are the 
variables.
}
\item{ina}{
A numerical vector with 1s, 2s, 3s and so one indicating the two groups. Be careful, the function is desinged to 
accept numbers greater than zero. 
Alternatively it can be a factor variable.
}
\item{type}{
This is for the variances test and can be either "levene" or "bf" corresponding to Levene's or Brown-Forsythe's 
testing procedure.
}
\item{treat}{
In the case of the blocking ANOVA this argument plays the role of the "ina" argument.
}
\item{block}{This item, in the blocking ANOVA denotes the subjects which are the same. 
Similarly to "ina" a numeric vector with 1s, 2s, 3s and so on.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
The Welch's F-test (without assuming equal variances) is performed with the "ftests" function. The "anovas" 
function perform the classical (Fisher's) one-way analysis of variance (ANOVA) which assumes equal variance 
across the groups. 

The "vartests" perform hypothesis test for the equality of the variances in two ways, either via the Levene or via 
the Brown-Forshythe procedure. Levene's test employs the means, whereas the Brown-Forsythe procedure employs the 
medians and is therefore more robust to outliers. The "var2tests" implement the classical F test.

The "block.anova" is the ANOVA with blocking, randomised complete block design (RCBD). In this case, for every 
combination of the block and treatment values, there is only one observation. The mathematics are the same as 
in the case of two way ANOVA, but the assumptions different and the testing procedure also different. 
In addition, no interaction is present. 
}

\value{
A matrix with the test statistic and the p-value of each test. 
}

\references{
B.L. Welch (1951). On the comparison of several mean values: an alternative approach. Biometrika, 38(3/4), 330-336.

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
\code{\link{ttests}
}
}

\examples{
x <- matrix( rnorm(300 * 50), ncol = 50 )
## 300 observations in total
ina <- rbinom(300, 3, 0.6) + 1   
a1 <- ftests(x, ina) 
a2 <- anovas(x, ina) 
a3 <- vartests(x, ina) 
x <- NULL
}

\keyword{ F-tests }
\keyword{ analysis of variance }
\keyword{ equality of variances }
\keyword{ blocking ANOVA }
