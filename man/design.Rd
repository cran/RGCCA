\name{design}
\alias{design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Design matrix}
\description{The utilisation of RGCCA requires the construction of the design matrix which characterize the connection between blocks. 
Elements of the (symmetric) design matrix C = (c_jk) is equal to 1 if block j and block k are connected and 0 otherwise.}
\usage{
design(A)
}
\arguments{\item{A}{a R-list containg the blocks to include in the RGCCA analysis}}

\value{\item{C}{Design matrix}}

\author{Arthur Tenenhaus <arthur.tenenhaus@supelec.fr>}

\examples{
if(interactive()){
data(Russett)
X_agric =as.matrix(Russett[ ,c("gini","farm","rent")])
X_ind = as.matrix(Russett[ ,c("gnpr","labo")])
X_polit = as.matrix(Russett[ ,c("demostab", "dictatur")])
A = list(X_agric, X_ind, X_polit)
design(A) # output: C
}
}
