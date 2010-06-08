\name{tau.estimate}
\alias{tau.estimate}
\title{Optimal Shrinkage intensity parameters}
\description{tau.estimate() estimate the shrinkage intensity parameter 
using the analytical formula proposed in Schafer and Strimmer (2005) for computing a 
regularized covariance matrix (using target A by default). It is possbile to specify another target.}
\usage{tau.estimate(X, tar = diag(ncol(X)))}
\arguments{\item{X}{Dataset on which estimated covariance matrix}
\item{tar}{target matrix, by defualt being eqal to the identity}
}
\value{\item{tau}{Optimal shrinkage intensity parameter}}
\references{Schafer J. and Strimmer K. (2005): A shrinkage approach to large-scale 
covariance matrix estimation and implications for functional genomics, Statistical 
applications in genetics and molecular biology 4, 1, Article 32.

Jelizarow M., Guillemot V., Tenenhaus A., Strimmer K., Boulesteix A.-L., (2010):
Over-optimism in bioinformatics: an illustration. Bioinformatics 26:1990-1998
}
\author{
Vincent Guillemot <vincent.guillemot@ibe.med.uni-muenchen.de>  and
 
Monika Jelizarow <jelizarow@ibe.med.uni-muenchen.de>  and           

Arthur Tenenhaus <arthur.tenenhaus@supelec.fr>
}

\examples{
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
A = list(X_agric, X_ind, X_polit)
tau = unlist(lapply(A, tau.estimate))
}
