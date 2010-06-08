rgcca <-
function(A, C, tau = "optimal", scheme = "centroid", scale = TRUE, layout = TRUE){

######################
# A.	INITIALISATION #
######################
#########################
## Data standardization #
#########################          
nbloc = length(A)

### standardised step
A = lapply(A, as.matrix)
if (scale == TRUE) {A = lapply(A, myscale)}
   
if (!is.numeric(tau)) {
print("Optimal Shrinkage intensity paramaters are estimated")
tau = unlist(lapply(A, tau.estimate))
}

else { 
if (is.numeric(tau)){print("Shrinkage intensity paramaters are chosen manually")}
}

#	A1.	Choose J arbitrary vectors \sim{a}_1^(0), \sim{a}_2^(0), ..., \sim{a}_J^(0) 
a = list()
	 for (q in 1:nbloc)
	   a[[q]] = rnorm(ncol(A[[q]]))

           
 
#	A2.	Compute vectors a_1^(0), a_2^(0), ..., a_J^(0) verifying the constraints of the general otpimization problem
#	and A3.	Compute the outer components
M = list()
Y = matrix(0, nrow(A[[1]]), nbloc)
for (q in 1:nbloc){
          M[[q]] = solve(tau[q]*diag(ncol(A[[q]])) + ((1-tau[q])/nrow(A[[q]]))*(t(A[[q]])%*%A[[q]]))
          a[[q]] = as.vector(1/sqrt(t(a[[q]])%*%M[[q]]%*%a[[q]]))*M[[q]]%*%a[[q]]
		      Y[ ,q] = A[[q]] %*% a[[q]]
}
	 
iter = 1
crit = numeric()
converg = numeric()
Z = matrix(0, nrow(A[[1]]), nbloc)
AVE_X = numeric()
a_temp = a

#-------------------------------------------------------

if ((scheme != "horst" ) & (scheme != "factorial") & (scheme != "centroid"))
{cat("ERROR : choose one of the three following schemes : horst, centroid or factorial")}
else{
if (scheme == "horst"){print("Computation of the RGCCA block components based on the Horst scheme")}
if(scheme == "factorial"){print("Computation of the RGCCA block components based on the Factorial scheme")}
if(scheme == "centroid"){print("Computation of the RGCCA block components based on the Centroid scheme")}
#-------------------------------------------------------
    repeat{
      Yold = Y

##########################################
# B.	INNER COMPONENT FOR THE BLOCK X_q  #
##########################################

#-------------------------------------------------------
    if (scheme == "horst"){
      for (q in 1:nbloc){
        #	Compute the inner component based on the Horst scheme:
        Z[,q] = rowSums(matrix(rep(C[q, ], nrow(A[[q]])), nrow(A[[q]]), nbloc, byrow = TRUE)*Y)

#########################################        
# C.	OUTER COMPONENT FOR THE BLOCK X_i #
#########################################   
        #	C.1.	Compute the outer weight
        a[[q]] = as.vector(1/sqrt(t(Z[, q])%*%A[[q]]%*%M[[q]]%*%t(A[[q]])%*%Z[, q]))*(M[[q]]%*%t(A[[q]])%*%Z[, q])
        #	C2.	Compute the outer component
        Y[, q] = A[[q]]%*%a[[q]]        
      }
    }

#-------------------------------------------------------
    if (scheme == "factorial"){
      for (q in 1:nbloc){
        #	Compute the inner component based on the factorial scheme:
        Z[, q] = rowSums(matrix(rep(C[q, ], nrow(A[[q]])), nrow(A[[q]]), nbloc, byrow = TRUE)*matrix(rep(mycov(Y[, q], Y), nrow(A[[q]])), nrow(A[[q]]), nbloc, byrow = TRUE)*Y)
        
#########################################        
# C.	OUTER COMPONENT FOR THE BLOCK X_i #
#########################################   
        #	C.1.	Compute the outer weight
        a[[q]] = as.vector(1/sqrt(t(Z[, q])%*%A[[q]]%*%M[[q]]%*%t(A[[q]])%*%Z[, q]))*(M[[q]]%*%t(A[[q]])%*%Z[, q])
        #	C2.	Compute the outer component
        Y[, q] = A[[q]]%*%a[[q]]
      }
    }

#-------------------------------------------------------
    if (scheme == "centroid"){
      for (q in 1:nbloc){
        #	Compute the inner component based on the centroid scheme:
        Z[,q] = rowSums(matrix(rep(C[q, ], nrow(A[[q]])), nrow(A[[q]]), nbloc, byrow = TRUE)*sign(matrix(rep(mycov(Y[, q], Y), nrow(A[[q]])), nrow(A[[q]]), nbloc, byrow = TRUE))*Y)
        
#########################################        
# C.	OUTER COMPONENT FOR THE BLOCK X_i #
#########################################   
        #	C.1.	Compute the outer weight
        a[[q]] = as.vector(1/sqrt(t(Z[, q])%*%A[[q]]%*%M[[q]]%*%t(A[[q]])%*%Z[, q]))*(M[[q]]%*%t(A[[q]])%*%Z[, q])
        #	C2.	Compute the outer component
        Y[, q] = A[[q]]%*%a[[q]]
      }
    }


    num_converg <- sum((rowSums(Yold)-rowSums(Y))^2)
    den_converg <- sum(rowSums(Yold)^2)
    converg[iter] <- num_converg/den_converg
      
    #check for convergence of the RGCCA alogrithm to a fixed point of the stationnary equations
    stationnary_point = rep(FALSE, length(A))
    for (ind in 1:length(A)){
        stationnary_point[ind] = sum(round(abs(a_temp[[ind]]-a[[ind]]), 8) < .Machine$double.eps) == ncol(A[[ind]])
    }
    a_temp = a


    if (scheme == "horst"){crit[iter] = sum(C*mycov(Y))}
    if(scheme == "factorial"){crit[iter] = sum((C*mycov(Y)^2))}
    if(scheme == "centroid"){crit[iter] = sum(C*abs(mycov(Y)))}
    
    if ((converg[iter] <.Machine$double.eps & sum(stationnary_point)== length(A)) | iter > 10000)  {break}
    iter = iter + 1
 }
 
     if (sum(stationnary_point) == length(A)){print("The RGCCA algorithm converged to a fixed point of the stationary equations")}
     else {print("The RGCCA algorithm did not converge after 10000 iterations.")}

  
    if (layout == TRUE){
    plot(1:length(which(crit != 0)), log(crit[which(crit != 0)]), xlab = "iteration", ylab = "criteria") 
    x11(); pairs(Y, labels = paste("Y", 1:ncol(Y)))
    }
    
    lambda = numeric(nbloc)
    for (q in 1:nbloc) {
      lambda[q] = mycov(Y[, q], Z[, q]) 
     
    # AVE for each block
      AVE_X[q] = mean(cor(A[[q]], Y[, q])^2)
    }
    # AVE outer model 
    AVEouter = sum((unlist(lapply(A, function(x) ncol(x)))* AVE_X))/sum(unlist(lapply(A, function(x) ncol(x))))
    # AVE inner model
    AVEinner =  sum(C*cor(Y)^2/2)/(sum(C)/2)
    AVE = list(AVE_per_block = AVE_X, AVE_outer_model = AVEouter, AVE_inner_model = AVEinner)
                         
    result = list(Y = Y, Z = Z, a = a, crit = crit[which(crit != 0)], converg = converg[which(converg != 0)], lambda = lambda, AVE = AVE, C = C, tau = tau, scheme = scheme)
    return(result)
     
}                                                                             
}

