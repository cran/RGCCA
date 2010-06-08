design <-
function(A){
C = matrix(0, length(A), length(A))
colnames(C) = rownames(C) = paste("X", 1:length(A), sep ="")
data.entry(C)
}

