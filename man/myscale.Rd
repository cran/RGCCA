\name{myscale}
\alias{myscale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standardization (zero mean and unit variance) of matrix-like objects }
\description{myscale is function whose default method centers and scales the columns of a numeric matrix.}
\usage{myscale(A, center = TRUE, scale = TRUE)}
\arguments{
  \item{A}{a numeric matrix }
  \item{center}{a logical value. If center = TRUE, each column is transformed to have zero mean}
  \item{scale}{a logical value. If scale = TRUE, each column is transformed to have unit variance}
}

\value{
\item{A}{The centered and/or scaled matrix. The numeric centering and scalings used (if any) are returned as attributes "scaled:center" and "scaled:scale"}
}

\author{Arthur Tenenhaus <arthur.tenenhaus@supelec.fr>}

\examples{
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
A = list(X_agric, X_ind, X_polit)
A = lapply(A, myscale)
}
