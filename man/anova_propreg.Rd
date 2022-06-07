\name{Significance testing for the coefficients of Quasi binomial or the quasi Poisson regression}
\alias{anova_propreg}
\alias{anova_qpois.reg}
\title{
Significance testing for the coefficients of Quasi binomial or the quasi Poisson regression
}

\description{
Significance testing for the coefficients of Quasi binomial or the quasi Poisson regression.
}

\usage{
anova_propreg(mod, poia = NULL)
anova_qpois.reg(mod, poia = NULL)
}

\arguments{
\item{mod}{
An object as returned by the "prop.reg" or the "qpois.reg" function.
}
\item{poia}{
If you want to test the significance of a single coefficient this must be a number. In this case, the "prop.reg" or the "qpois.reg" function contains
this information. If you want more coefficients to be testes simultaneously, e.g. for a categorical predictor, then this must 
contain the positions of the coefficients. If you want to see if all coefficients are zero, like an overall F-test, leave this NULL. 
}
}

\details{
Even though the name of this function starts with anova it is not an ANOVA type significance testing, but a Wald type.  
}

\value{
A vector with three elements, the test statistic value, its associated p-value and the relevant degrees of freedom.
}

\references{
Papke L. E. & Wooldridge J. (1996). Econometric methods for fractional response variables with 
an application to 401(K) plan participation rates. Journal of Applied Econometrics, 11(6): 619-632.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{prop.reg}, \link{qpois.reg}, \link{univglms}, \link{score.glms}, \link{logistic_only}
}
}

\examples{
\donttest{
y <- rbeta(1000, 1, 4)
x <- matrix(rnorm(1000 * 3), ncol = 3)
a <- prop.reg(y, x)
## all coefficients are tested
res<-anova_propreg(a)  
## the first predictor variable is tested 
res<-anova_propreg(a, 2)  
a  ## this information is already included in the model output
## the first and the second predictor variables are tested 
res<-anova_propreg(a, 2:3)
}
}

\keyword{ Significance testing }
\keyword{ GLMs }
\keyword{ Quasi regression }
\keyword{ fractional response }