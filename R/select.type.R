select.type <- function(type, tau, C, scheme, A, ncomp){
  
  # This function translates the type string into the appropriate tau and function
  # Predefined functions:
  # c("rgcca",
  #  "cpca-w", "gcca", "hpca", "maxbet-b", "maxbet", "maxdiff-b", "maxdiff",
  #  "maxvar-a", "maxvar-b", "maxvar", "niles", "r-maxvar", "rcon-pca",
  #  "ridge-gca", "sabscor", "ssqcor", "ssqcor", "ssqcov-1", "ssqcov-2",
  #  "ssqcov", "sum-pca", "sumcor", "sumcov-1", "sumcov-2", "sumcov")
  
  J = length(A)
  C.hierarchical <- matrix(0, J+1, J+1)
  C.hierarchical[1:J, J+1] = C.hierarchical[J+1, 1:J] = 1
 

  if (tolower(type) == "rgcca"){
    scheme   <- scheme
    tau      <- tau
    C        <- C
    A        <- A
    ncomp    <- ncomp
  } 
  
  else if (tolower(type) == "sumcor"){
    scheme   <- function(x) x
    tau      <- rep(0, J)
    C        <- 1-diag(J)
    A        <- A
    ncomp    <- ncomp
  } 
  
  else if (tolower(type) == "ssqcor"){
    scheme   <- function(x) x^2
    tau      <- rep(0, J)
    C        <- 1-diag(J)
    A        <- A
    ncomp    <- ncomp
  } 
  
  else if (tolower(type) == "sabscor"){
    scheme   <- "centroid"
    tau      <- rep(0, J)
    C        <- 1-diag(J)
    A        <- A
    ncomp    <- ncomp
  } 
  
  else if (tolower(type)%in%c("sumcov","sumcov-1", "maxbet")){
    scheme   <- function(x) x
    tau      <- rep(1, J)
    C        <- matrix(1, J, J)
    A        <- A
    ncomp    <- ncomp
  } 
  
  else if (tolower(type)%in%c("sumcov-2", "maxdiff")){
    scheme   <- function(x) x
    tau      <- rep(1, J)
    C        <- 1 - diag(J)
    A        <- A
    ncomp    <- ncomp
  } 
  
  else if (tolower(type)%in%c("ssqcov","ssqcov-1", "maxbet-b")){
    scheme   <- function(x) x^2
    tau      <- rep(1, J)
    C        <- matrix(1, J, J)
    A        <- A
    ncomp    <- ncomp
  } 
  
  else if (tolower(type)%in%c("ssqcov-2", "maxdiff-b")){
    scheme   <- function(x) x^2
    tau      <- rep(1, J)
    C        <- 1 - diag(J)
    A        <- A
    ncomp    <- ncomp
  } 
  
  else if (tolower(type) == "rcon-pca"){
    scheme   <- scheme
    tau      <- tau
    C        <- C.hierarchical
    A[[J+1]] <- Reduce("cbind", A)
    ncomp    <- c(ncomp, ncomp[1])
  } 
  
  else if (tolower(type)%in%c("maxvar-b", "gcca", "niles", "maxvar")){
    scheme   <- function(x) x^2
    tau      <- rep(0, J+1)
    C        <- C.hierarchical
    A[[J+1]] <- Reduce("cbind", A)
    ncomp    <- c(ncomp, ncomp[1])
  } 
  
  else if (tolower(type)%in%c("maxvar-a", "cpca-w", "sum-pca")){
    scheme   <- function(x) x^2
    tau      <- c(rep(1, J), 0)
    C        <- C.hierarchical
    A[[J+1]] <- Reduce("cbind", A)
    ncomp    <- c(ncomp, ncomp[1])
  } 
  
  else if (tolower(type) == "ridge-gca"){
    scheme   <- function(x) x^2
    tau      <- c(tau[1:J], 0)
    C        <- C.hierarchical
    A[[J+1]] <- Reduce("cbind", A)
    ncomp    <- c(ncomp, ncomp[1])
  } 
  
  else if (tolower(type) == "r-maxvar"){
    scheme   <- function(x) x^2
    tau      <- tau
    C        <- C.hierarchical
    A[[J+1]] <- Reduce("cbind", A)
    ncomp    <- c(ncomp, ncomp[1])
  } 
  
  else if (tolower(type) == "hpca"){
    scheme   <- function(x) x^4
    tau      <- c(rep(1, J), 0)
    C        <- C.hierarchical
    A[[J+1]] <- Reduce("cbind", A)
    ncomp    <- c(ncomp, ncomp[1])
  }
  
  return(list(scheme = scheme, tau = tau, C = C, A = A))
}
