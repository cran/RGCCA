\name{RGCCA-package}
\alias{RGCCA-package}
\alias{RGCCA}
\docType{package}
\title{RGCCA Package}
\description{ This package supplies an efficient methodology to unravel relationships between several blocks 
of variables observed on the same set of individuals. Regularized Generalized Canonical Correlation Analysis (RGCCA) 
is a generalization of regularized canonical correlation analysis to three or more sets of 
variables and constitutes a general framework for many multi-block data analysis methods.  
It combines the power of multi-block data analysis methods (maximization of well identified criteria) 
and the flexibility of PLS path modeling (the researcher decides which blocks are connected and which are not).  
Searching for a fixed point of the stationary equations related to the RGCCA optimization problem, a new monotone convergent 
algorithm is proposed in this package.
}
\details{
\tabular{ll}{
Package: \tab RGCCA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-06-08\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}
\author{Arthur Tenenhaus - Maintainer: Arthur Tenenhaus <Arthur.Tenenhaus@supelec.fr>}
\references{Tenenhaus A. and Tenenhaus M., Regularized Generalized Canoncial Correlation Anlaysis, Psychometrika, 2010}
\keyword{package}
\examples{
#load Russett data
data(Russett)
# Construction of the three blocks of variables observed on the same set 
# of individuals. These blocks are stored in a R list.
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
A = list(X_agric, X_ind, X_polit)

#Define the design matrix (output = C) 
if(interactive()){design(A)}
if(interactive() == FALSE){C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)}

# For each block, the optimal shrinkage constant parameter is determined 
# using the Schafer and Strimmer's formula (tau = "optimal"). 
result.rgcca = rgcca(A, C, tau = "optimal", "factorial", scale = TRUE, 
                     layout = TRUE)  

# Instead of choosing to compute the optimal shrinkage intensity 
# parameters, we can define them manually:
result.rgcca = rgcca(A, C, tau = c(1, 1, 0), "factorial", 
                     scale = TRUE, layout = TRUE)

lab = as.vector(apply(Russett[, 10:12], 1, which.max))
x11();plot(result.rgcca$Y, col = "white", xlab = "Y1 (Agric. inequality)", 
           ylab = "Y2 (Industrial Development)")
text(result.rgcca$Y[, 1], result.rgcca$Y[, 2], Russett[, 1], col = lab)
}

