#' The function sgccak() is called by sgcca() and does not have to be used by
#' the user. sgccak() enables the computation of SGCCA block components, outer
#' weight vectors, etc., for each block and each deflation stage.
#' @inheritParams sgcca
#' @inheritParams rgccak
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{crit}{A vector of integer that contains for each component
#' the values of the analysis criteria across iterations.}
#' @title Internal function for computing the SGCCA parameters (SGCCA block
#' components, outer weight vectors etc.)
#' @noRd
sgccak <- function(A, C, sparsity = rep(1, length(A)),
                   scheme = "centroid", tol = 1e-08,
                   init = "svd", bias = TRUE, verbose = FALSE,
                   na.rm = TRUE, response = NULL,
                   disjunction = FALSE, n_iter_max = 1000) {
  if (is.function(scheme)) {
    g <- scheme
  } else {
    switch(scheme,
      "horst" = {
        g <- function(x) x
      },
      "factorial" = {
        g <- function(x) x^2
      },
      "centroid" = {
        g <- function(x) abs(x)
      }
    )
  }

  dg <- Deriv::Deriv(g, env = parent.frame())

  ### Initialization
  init_object <- sgcca_init(
    A, init, bias, na.rm, sparsity, response, disjunction
  )
  a <- init_object$a
  Y <- init_object$Y

  iter <- 1
  crit <- NULL
  crit_old <- sum(C * g(cov2(Y, bias = bias)))
  a_old <- a

  repeat {
    update_object <- sgcca_update(
      A, bias, na.rm, sparsity, response, disjunction, dg, C, a, Y, init_object
    )
    a <- update_object$a
    Y <- update_object$Y

    # Print out intermediate fit
    crit <- c(crit, sum(C * g(cov2(Y, bias = bias))))

    if (verbose) {
      cat(
        " Iter: ", formatC(iter, width = 3, format = "d"),
        " Fit: ", formatC(crit[iter], digits = 8, width = 10, format = "f"),
        " Dif: ", formatC(crit[iter] - crit_old,
          digits = 8, width = 10, format = "f"
        ), "\n"
      )
    }
    stopping_criteria <- c(
      drop(crossprod(unlist(a, FALSE, FALSE) - unlist(a_old, FALSE, FALSE))),
      abs(crit[iter] - crit_old)
    )

    if (any(stopping_criteria < tol) || (iter > n_iter_max)) {
      break
    }

    crit_old <- crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  if (iter > n_iter_max) {
    stop_rgcca(
      "The SGCCA algorithm did not converge after ", n_iter_max,
      " iterations."
    )
  }
  if (verbose) {
    if (iter <= n_iter_max) {
      message(
        "The SGCCA algorithm converged to a stationary point after ",
        iter - 1, " iterations \n"
      )
    }
    plot(crit, xlab = "iteration", ylab = "criteria")
  }

  result <- sgcca_postprocess(
    A, a, Y, g, na.rm, sparsity, tol, response, disjunction
  )
  return(list(Y = result$Y, a = result$a, crit = crit))
}
