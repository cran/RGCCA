mycov <-
function(x, y = NULL){
if (is.null(y)){
  if (is.vector(x)){n = length(x);((n-1)/n)*var(x)}
  else{n = nrow(x);((n-1)/n)*cov(x)}
}
else {((length(x)-1)/length(x))*cov(x, y)}
}

