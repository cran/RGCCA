\name{Russett}
\alias{Russett}
\docType{data}
\title{Russett data}
\description{
The Russett data set (Russett, 1964) are studied in Gifi (1990). Three blocks of variables have been defined for 47 countries. The first block X1 = [GINI, FARM, RENT] is related to "Agricultural Inequality". The second block X2 = [GNPR, LABO] describes "Industrial Development". The third one X3 = [INST, ECKS, DEAT] measures "Political Instability". An additional variable DEMO describes the political regime: stable democracy, unstable democracy or dictatorship. Russett collected this data to study relationships between Agricultural Inequality, Industrial Development and Political Instability. Russett's hypotheses can be formulated as follows: It is difficult for a country to escape dictatorship when its agricultural inequality is above-average and its industrial development below-average.
}
\usage{data(Russett)}
\format{
  A data frame with 47 observations on the following 11 numeric variables.
  \describe{
    \item{\code{gini}}{Inequality of land distribution}
    \item{\code{farm}}{\% farmers that own half of the land }
    \item{\code{rent}}{\% farmers that rent all their land}
    \item{\code{gnpr}}{Gross national product per capita ($1955)}
    \item{\code{labo}}{\% of labor force employed in agriculture}
    \item{\code{inst}}{Instability of executive (45-61)}
    \item{\code{ecks}}{Number of violent internal war incidents (46-61)}
    \item{\code{death}}{Number of people killed as a result of civic group violence (50-62)}
    \item{\code{demostab}}{binary variable equal to 1 for stable democraty and 0 otherwise}
    \item{\code{demoinst}}{binary variable equal to 1 for unstable democraty and 0 otherwise}
    \item{\code{dictator}}{binary variable equal to 1 for dictatorship and 0 otherwise}
  }
}
\references{Russett B.M. (1964), Inequality and Instability: The Relation of Land Tenure to Politics, World Politics 16:3, 442-454.

Gifi, A. (1990), Nonlinear multivariate analysis, Chichester: Wiley.
}
\examples{
#Loading of the Russett dataset
data(Russett)
#Russett is partitioned into three blocks (X_agric, X_ind, X_polit)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("inst", "ecks", "death", "demostab", 
                                 "demoinst", "dictator")])
A = list(X_agric, X_ind, X_polit)
lapply(A, dim)
}
\keyword{datasets}
