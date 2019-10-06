\name{MLE of the ordinal model without covariates}
\alias{ordinal.mle}
\title{
MLE of the ordinal model without covariates
}

\description{
MLE of the ordinal model without covariates.
}

\usage{
ordinal.mle(y, link = "logit")
}

\arguments{
\item{y}{
A numerical vector with values 1, 2, 3,..., not zeros, or an ordered factor. 
}
\item{link}{
This can either be "logit" or "probit". It is the link function to be used. 
}
}

\details{
Maximum likelihood of the ordinal model (proportional odds) is implemented. See for example the "polr" command 
in R or the examples.
}

\value{
A list including:
\item{loglik}{
The log-likelihood of the model.
}
\item{a}{
The intercepts (threshold coefficients) of the model.
}
}

\references{
Agresti, A. (2002) Categorical Data. Second edition. Wiley. 
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{beta.mle}, \link{diri.nr2}
}
}

\examples{
y <- factor( rbinom(100,3,0.5), ordered = TRUE )
res<-ordinal.mle(y)
res<-ordinal.mle(y, link = "probit")
}

\keyword{ Ordinal model }
\keyword{ proportional odds }

