#' @export
#' @rdname print
#' @order 5
print.rgcca_stability <- function(x, ...) {
  stopifnot(is(x, "rgcca_stability"))
  cat("RGCCA stability object obtained with", x$n_boot, "bootstrap samples.")
  cat("\n")
}
