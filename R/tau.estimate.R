tau.estimate <- function(X,tar = diag(ncol(X))) {
# check for possibly incorrect inputs
    if (is.matrix(X)==TRUE && is.numeric(X)==FALSE) {                           
    stop("The data matrix must be numeric!")}
    # correlation matrix
    p = ncol(X) ; n = nrow(X)
    corm = cor(X)                                                               
    # Standardization
    X = myscale(X)
    # matrix of the var_hat_s[i,j] in the formula for lambda                                                             
    v = matrix(0,ncol=p,nrow=p)                                                 
                                                                                
    ind = which(upper.tri(v),arr.ind=TRUE)
    vtmp = apply(ind,1, function(ii) sum( ( X[,ii[1]]*X[,ii[2]]- mean(X[,ii[1]]*X[,ii[2]]) )**2 ) )
    v[ind] = (n/((n-1)^3))*vtmp ; 
    v = v + t(v)    
    
    # matrix of the f_[ij]'s in the formula for lambda
    f = matrix(0,ncol=p,nrow=p)                                                 
                                                                                
    ind = which(upper.tri(f) & tar != 0,arr.ind=TRUE)
    ftmp = apply(ind,1, function(ii) 
                  1/2*sum((X[,ii[1]]*X[,ii[1]]- mean(X[,ii[1]]*X[,ii[1]]) 
                  + X[,ii[2]]*X[,ii[2]]- mean(X[,ii[2]]*X[,ii[2]]))
                  * (X[,ii[1]]*X[,ii[2]]- mean(X[,ii[1]]*X[,ii[2]])))
                )
    f[ind] = (n/((n-1)^3))*ftmp ; f = f + t(f)    

    corapn =  cov2cor(tar)
    # matrix of the single values of the denominator in the formula for lambda
    d = (corm - corapn)^2                                                       
                                                                             
    # formula for the shrinkage intensity (see Schäfer&Strimmer (2005))
    lambda = (sum(v)- sum(corapn*f))/sum(d)                                     
    lambda = max(min(lambda, 1), 0)                                              

    return(lambda)


 }