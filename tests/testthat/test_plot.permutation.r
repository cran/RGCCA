#' # test plot.permutation
#-------------------------
set.seed(0)
data(Russett)
A <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
perm.out <- rgcca_permutation(A,
  par_type = "tau", n_perms = 5,
  n_cores = 1, verbose = FALSE
)

A2 <- c(A, A)
names(A2) <- NULL
perm.out2 <- rgcca_permutation(A2,
  par_type = "sparsity", n_perms = 2,
  n_cores = 1, par_length = 4, verbose = FALSE
)

test_that("plot.permutation produces the expected crit plot", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "Permutation crit", plot.permutation(perm.out, type = "crit")
  )
})

test_that("plot.permutation produces the expected zstat plot", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "Permutation zstat", plot.permutation(perm.out, type = "zstat")
  )
})

test_that("plot.permutation produces the expected zstat plot with legend", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "Permutation legend", plot.permutation(
      perm.out,
      type = "zstat", show_legend = TRUE
    )
  )
})

test_that("plot.permutation produces the expected crit plot with many blocks", {
  skip_on_cran()
  vdiffr::expect_doppelganger(
    "Permutation many blocks", plot.permutation(perm.out2, type = "crit")
  )
})
