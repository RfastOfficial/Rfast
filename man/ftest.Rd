\name{Multi-sample tests for vectors}
\alias{ftest}
\alias{anova1}
\alias{kruskaltest}
\alias{var2test}
\alias{ttest2}
\alias{mcnemar}
\alias{cqtest}
\alias{block.anova}
\alias{twoway.anova}
\title{
Multi-sample tests for vectors
}

\description{
Multi-sample tests for vectors.
}

\usage{
ftest(x, ina, logged = FALSE)
anova1(x, ina, logged = FALSE)
kruskaltest(x, ina, logged = FALSE) 
var2test(x, y, alternative = "unequal", logged = FALSE)
mcnemar(x, y, logged = FALSE)
ttest2(x, y, paired = FALSE, logged = FALSE)
cqtest(x, treat, block, logged = FALSE) 
block.anova(x, treat, block, logged = FALSE) 
twoway.anova(y, x1, x2, interact = FALSE, logged = FALSE)
}

\arguments{
\item{x}{
A numerical vector with the data.
}
\item{y}{
A numerical vector with the data.
}
\item{ina}{
A numerical vector with 1s, 2s, 3s and so one indicating the two groups. Be careful, the function is desinged to 
accept numbers greater than zero. 
Alternatively it can be a factor variable.
}
\item{paired}{
This is for the two sample t-test  only and is TRUE or FALSE specifying whether the two samples are paired or not.
}
\item{alternative}{
This can either be "unequal", "greater" or "less".
}
\item{treat}{
In the case of the blocking ANOVA and Cochran's Q test, this argument plays the role of the "ina" argument.
}
\item{block}{This item (in the blocking ANOVA and Cochran's Q test) denotes the subjects which are the same. 
Similarly to "ina" a numeric vector with 1s, 2s, 3s and so on.
}
\item{x1}{
The first factor in the two way ANOVA.
}
\item{x2}{
The second factor in the two way ANOVA. The orderis not important.
}
\item{interact}{
Should interaction in the two way ANOVA be included? The default value is FALSE (no interaction). 
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
The Welch's F-test (without assuming equal variances) is performed with the "ftest" function. The "anova" function 
perform the classical (Fisher's) one-way analysis of variance (ANOVA) which assumes equal variance across the groups. 
The "kruskaltest" performs the Kruskal-Wallis non parametric alternative to analysis of variance test.
The "var2tests" implement the classical F test for the equality of two sample variances. The "cqtest" performs the 
Cocrhan's Q test for the equality of more than two groups whose values are strictly binary (0 or 1). This is a 
generalisation of the McNemar's test in the multi-sample case. The "block.anova" is the ANOVA with blocking, 
randomised complete block design (RCBD). In this case, for every combination of the block and treatment values, there 
is only one observation. The mathematics are the same as in the case of "twoway.anova", but the assumptions different
and the testing procedure also different. In addition, no interaction is present. 
}

\value{
A vector with the test statistic and the p-value of each test. For the case of t-test, an extra column with the 
degrees of freedom is given. For the two way ANOVA there can can be either 2 or three F test statistics and hence 
the same number of p-values. 
}

\references{
B.L. Welch (1951). On the comparison of several mean values: an alternative approach. Biometrika, 38(3/4), 330-336.

D.C. Montgomery (2001). Design and analysis of experiments (5th Edition). New York: John Wiley \& Sons.

McNemar Q. (1947). Note on the sampling error of the difference between correlated proportions or percentages. 
Psychometrika. 12(2):153-157. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ttests}, \link{ftests}
}
}

\examples{
x <- rnorm(200)
ina <- rbinom(200, 3, 0.5) + 1
res<-anova1(x, ina)
res<-ftest(x, ina)
ina <- rbinom(200, 1, 0.5) + 1
x1 <- x[ ina == 1 ]   ; x2 <- x[ ina == 2 ]
res<-ttest2(x1, x2)
res<-var2test(x1, x2)

## RCBD example 4.1 from Montgomery (2001), page 131-132
x <- c(9.3, 9.4, 9.2, 9.7, 9.4, 9.3, 9.4, 9.6, 9.6, 9.8, 9.5, 10,
10, 9.9, 9.7, 10.2)
tr <- rep(1:4, 4)
bl <- rep(1:4, each = 4)
res<-block.anova(x, tr, bl)
}

\keyword{ F-test }
\keyword{ analysis of variance }
\keyword{ equality of variances }
\keyword{ James test }
\keyword{ dependent binary data }
\keyword{ blocking ANOVA }
