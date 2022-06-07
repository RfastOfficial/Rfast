\name{Analysis of variance with a count variable}
\alias{poisson.anova}
\alias{geom.anova}
\alias{quasipoisson.anova}
\title{
Analysis of variance with a count variable
}

\description{
Analysis of variance with a count variable.
}

\usage{
poisson.anova(y, ina, logged = FALSE)
geom.anova(y, ina, type = 1, logged = FALSE)
quasipoisson.anova(y, ina, logged = FALSE)
}

\arguments{
\item{y}{
A numerical vector with discrete valued data, i.e. counts.
}
\item{ina}{
A numerical vector with discrete numbers starting from 1, i.e. 1, 2, 3, 4,... or a factor variable. 
This is suppose to be a categorical predictor. If you supply a continuous valued vector the function 
will obviously provide wrong results.
}
\item{type}{
This argument is for the geometric distribution. Type 1 refers to the case where the minimum is zero and type 2 for the case of the minimum being 1. 
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
This is the analysis of variance with Poisson or geometric distributed data. What we do is a log-likelihood ratio 
test. However, this is exactly the same as Poisson regression with a single predictor variable who happens to be 
categorical. Needless to say that this is faster function than the glm command in R. For the same purpose with 
a Bernoulli variable use \code{\link{g2Test}}. The quasinpoisson.anova is when in the glm function you specify 
family = quasipoisson. This is suitable for the case of over or under-dispersed data. 
}

\value{
A vector with two values, the difference in the deviances (or the scale difference in the case of quasi poisson)
and the relevant p-value. The quasipoisson.anova also returns the estimate of the \eqn{\phi} parameter.  
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{logistic.cat1}, \link{g2Test}, \link{poisson.anovas}, \link{anova}, \link{poisson_only}, 
\link{poisson.mle}
}
}

\examples{
y <- rpois(300, 10)
ina <- rbinom(300, 3, 0.5) + 1 
a1 <- poisson.anova(y, ina) 
a2 <- glm(y ~ ina, poisson) 

\donttest{
res<-anova(a2, test = "Chisq")
}
y <- rgeom(300, 0.7)
res<-geom.anova(y, ina) 
}

\keyword{ Analysis of variance }
\keyword{ Poisson distribution }
\keyword{ geometric distribution }
