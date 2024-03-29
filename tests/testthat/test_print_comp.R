data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

test_that("print_comp prints the expected text", {
  local_edition(3)
  expect_snapshot({
    fit.rgcca <- rgcca(blocks)
    res <- print_comp(fit.rgcca, n = 1, i = 1, outer = FALSE)
    cat(res)
  })
})

test_that("print_comp prints the expected text 2", {
  local_edition(3)
  expect_snapshot({
    fit.rgcca <- rgcca(blocks, ncomp = 2, sparsity = c(1, 1, 0.5))
    res <- print_comp(fit.rgcca, n = 2, i = 3, outer = FALSE)
    cat(res)
  })
})

test_that("print_comp prints the expected text 3", {
  local_edition(3)
  expect_snapshot({
    fit.rgcca <- rgcca(blocks)
    res <- print_comp(fit.rgcca, outer = TRUE)
    cat(res)
  })
})

test_that("print_comp prints the expected text 4", {
  local_edition(3)
  expect_snapshot({
    fit.rgcca <- rgcca(blocks, ncomp = 4, superblock = TRUE)
    res <- print_comp(fit.rgcca, outer = TRUE)
    cat(res)
  })
})
