myscale <-
function(A, center = TRUE, scale = TRUE){

if (center == TRUE & scale == TRUE){
# Centering step
A = scale(A, center = TRUE, scale = FALSE) 

# standardised step
std = sqrt(apply(A, 2, mycov))
A = A/matrix(rep(std, nrow(A)), nrow(A), ncol(A), byrow = TRUE)
attr(A, "std") = std
return(A)
}

if( center == TRUE & scale == FALSE) {A = scale(A, center = TRUE, scale = FALSE) ; return(A)} 

if( center == FALSE & scale == TRUE) {
std = sqrt(apply(A, 2, mycov))
A = A/matrix(rep(std, nrow(A)), nrow(A), ncol(A), byrow = TRUE)
attr(A, "std") = std
return(A)
} 
}

