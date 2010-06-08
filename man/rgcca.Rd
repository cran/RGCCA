\name{rgcca}
\alias{rgcca}
\title{Regularized Generalized Canonical Correlation Analysis (RGCCA)}
\description{Regularized Generalized Canonical Correlation Analysis (RGCCA) is a generalization of regularized 
canonical correlation analysis to three or more sets of variables. Given J matrices X1, X2, ..., XJ, which 
represent J sets of features observed on the same set of individuals (=> X1, ..., XJ must have same number of 
rows, but may (and usually will) have different numbers of columns), RGCCA aims at studying the relationships 
between these J blocks of variables. It constitutes a general framework for many multi-block data analysis methods.  
It combines the power of multi-block data analysis methods (maximization of well identified criteria) and the 
flexibility of PLS path modeling (the researcher decides which blocks are connected and which are not). Searching 
for a fixed point of the stationary equations related to the RGCCA optimization problem, this function implements 
a new monotone convergent algorithm (i.e. the bounded criteria to be maximized increases at each step of the iterative 
procedure) which is very similar to the PLS algorithm proposed by Herman Wold.  
}
\usage{rgcca(A, C, tau = "optimal", scheme = "centroid", scale = TRUE, 
     layout = TRUE)}
\arguments{
  \item{A}{a R-list containing the J blocks of variables X1, X2, ..., XJ}
  \item{C}{Design matrix describing relationships between blocks}
  \item{tau}{a 1*J vector containing the value for the shrinkage constant parameters tau_j, j=1, ... J. If tau = "optimal" (default value) the shrinkage intensity paramaters are estimated using the Schafer and Strimmer (2005) analytical formula}
  \item{scheme}{"horst", "factorial", "centroid"}
  \item{scale}{if scale = TRUE, each block are standardized to zero means and unit variances}
  \item{layout}{plot figures: (i) yaxis = Criterion to be maximized as a function of xaxis = the number of iteration. (ii) A matrix of scatterplots is produced for the matrix of outer components.}
}
\value{
     \item{Y}{a n*J matrix of RGCCA outer components}
     \item{Z}{a n*J matrix of RGCCA inner components}
     \item{a}{a R-list of outer weight vector}
     \item{crit}{value of the function to be optimized at the end of the procedure}
     \item{converg}{Rate of convergence of the algorithm}
     \item{lambda}{value for lambdas}
     \item{AVE}{Indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model)}
     \item{C}{Design matrix describing the relation between blocks}
     \item{tau}{a 1*J vector containing the value for the shrinkage intensity parameters tau_j, j=1, ... J}
     \item{scheme}{The scheme chosen by the user.}

}
\references{Tenenhaus A. and Tenenhaus M., Regularized Generalized Canoncial Correlation Anlaysis, Psychometrika, 2010}
\author{Arthur Tenenhaus <arthur.tenenhaus@supelec.fr>}
\seealso{tau.estimate}
\examples{
#############
# Example 1 #
#############
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
A = list(X_agric, X_ind, X_polit)

#Define the design matrix (output = C) 
if(interactive()){design(A)}
if(interactive() == FALSE){C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)}

result.rgcca = rgcca(A, C, tau = c(1, 1, 1), "factorial", scale = TRUE)

lab = as.vector(apply(Russett[, 10:12], 1, which.max))
x11();plot(result.rgcca$Y, col = "white", xlab = "Y1 (Agric. inequality)", 
           ylab = "Y2 (Industrial Development)")
text(result.rgcca$Y[, 1], result.rgcca$Y[, 2], Russett[, 1], col = lab)
 
#######################################
# example 2 : RGCCA and leave one out #
#######################################

Ytest = matrix(0, 47, 3)
tau = c(1, 1, 1)
for (i in 1:nrow(Russett)){
  B = lapply(A, function(x) x[-i, ])
  B = lapply(B, myscale)
  resB = rgcca(B, C, tau = tau, "factorial", scale = FALSE, layout = FALSE)
#  look for potential conflicting sign among components within the loo loop.
  for (k in 1:length(B)){
  if (cor(result.rgcca$a[[k]], resB$a[[k]]) >= 0) 
      resB$a[[k]] = resB$a[[k]] else resB$a[[k]] = -resB$a[[k]]
  }
  Btest =lapply(A, function(x) x[i, ])
  Btest[[1]]=(Btest[[1]]-attr(B[[1]],"scaled:center"))/(attr(B[[1]],"std"))
  Btest[[2]]=(Btest[[2]]-attr(B[[2]],"scaled:center"))/(attr(B[[2]],"std"))
  Btest[[3]]=(Btest[[3]]-attr(B[[3]],"scaled:center"))/(attr(B[[3]],"std"))

  Ytest[i, 1] = Btest[[1]]\%*\%resB$a[[1]]
  Ytest[i, 2] = Btest[[2]]\%*\%resB$a[[2]]
  Ytest[i, 3] = Btest[[3]]\%*\%resB$a[[3]]
}

lab = apply(Russett[, 10:12], 1, which.max)
plot(result.rgcca$Y, col = "white", xlab = "Y1 (Agric. inequality)", 
     ylab = "Y2 (Ind. Development)")
text(result.rgcca$Y[, 1], result.rgcca$Y[, 2], Russett[, 1], col = lab)
text(Ytest[, 1], Ytest[, 2], substr(Russett[, 1], 1, 1), col = lab)

}
