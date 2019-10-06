\name{ANOVA for two quasi Poisson regression models}
\alias{anova_quasipois.reg}
\title{
ANOVA for two quasi Poisson regression models
}

\description{
ANOVA for two quasi Poisson regression models.
}

\usage{
anova_quasipois.reg(mod0, mod1, n) 
}

\arguments{
\item{mod0}{
An object as returned by the "qpois.reg" function. This is the null model.
}
\item{mod1}{
An object as returned by the "qpois.reg" function. This is the alternative model.
}
\item{n}{
The sample size. This is necessary to calculate the degrees of freedom.
}
}

\details{
This is an ANOVA type significance testing for two quasi Poisson models. 
}

\value{
A vector with 4 elements, the test statistic value, its associated p-value and the relevant 
degrees of freedom of the numerator and the denominator.
}

\references{
Papke L. E. & Wooldridge J. (1996). Econometric methods for fractional response variables with 
an application to 401(K) plan participation rates. Journal of Applied Econometrics, 11(6): 619--632.

McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{anova_qpois.reg}, \link{qpois.reg}, \link{univglms}, \link{quasipoisson.anova} 
}
}

\examples{
y <- rnbinom(200, 10, 0.5)
x <- matrix(rnorm(200 * 3), ncol = 3)
a1 <- qpois.reg(x, y)
a0 <- qpois.reg(x[, 1], y)
res<-anova_quasipois.reg(a0, a1, 200)
b1 <- glm(y ~ x, family = quasipoisson)
b0 <- glm(y ~ x[, 1], family = quasipoisson)
res<-anova(b0, b1, test = "F")
c1 <- glm(y ~ x, family = poisson)
c0 <- glm(y ~ x[, 1], family = poisson)
res<-anova(c0, c1, test = "Chisq")
}
