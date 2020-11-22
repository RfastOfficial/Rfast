\name{Topological sort of a DAG}
\alias{topological_sort}
\title{
Topological sort of a DAG
}

\description{
Topological sort of a DAG.
}

\usage{
topological_sort(dag)
}

\arguments{
\item{dag}{
A square matrix representing a directed graph which contains 0s and 1s. If G[i, j] = 1 it 
means there is an arrow from node i to node j. When there is no edge between nodes i and j if G[i, j] = 0.  
}
}

\details{
The function is an R translation from an old matlab code.
}

\value{
A vector with numbers indicating the sorting. If the dag is not a Directed acyclic Graph, NA will be returned.
}

\references{
Chickering, D.M. (1995). A transformational characterization of equivalent Bayesian network structures. 
Proceedings of the 11th Conference on Uncertainty in Artificial Intelligence, Montreal, Canada, 87-98. 
}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{floyd}, \link{pc.skel} 
}
}

\examples{
G <- matrix(0, 5, 5)
G[2, 1] <- 1
G[3, 1] <- 1
G[4, 2] <- 1
G[5, 4] <- 1
res<-topological_sort(G)
G[2, 4] <- 1
res<-topological_sort(G)
}

\keyword{ Topological sort }
\keyword{ DAG }