\name{mycov}
\alias{mycov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standardization}
\description{
 mycov compute the variance of x and the covariance of x and y if these are vectors
}
\usage{
mycov(x, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, matrix or data frame}
  \item{y}{a numeric vector, matrix or data frame}
}
\details{
The denominator n is used (instead of n-1) which gives a biased estimator of the (co)variance for i.i.d. observations. 
}

\value{
\item{var/cov}{estimation of the variance (resp. covariance) of x (resp. x and y)}
}

\author{Arthur Tenenhaus <arthur.tenenhaus@supelec.fr>}



\examples{
x = rnorm(10)
var.x = mycov(x)
y = rnorm(10)
cov.xy = mycov(x, y)
}
