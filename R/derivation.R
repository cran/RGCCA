#' The function derivation() allows computing the derivative of the g function. 
#' The g function must be convex.
#' @param g a convex function 
#' @return \item{dg}{the derivate of the g function}
#' @title compute the derivative of the g function
#' @export derivation
#' @importFrom methods formalArgs
#' @importFrom utils capture.output
#' @importFrom stats deriv as.formula

derivation <- function(g){
  arg <- formalArgs(g)
  gchar <- gsub(sprintf("function\\(%s\\)", arg),"", capture.output(g))
  dg <- deriv(as.formula( paste("~", gchar[1]) ), arg) 
}