# '# Test intersection_list
#
# '''
# Load data
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11],
  target = Russett[, 11]
)

# Add missing values
blocks[[1]][c(2, 4, 8), ] <- NA
blocks[[2]][c(12, 23), 1] <- NA
blocks[[2]][17, 2] <- NA
blocks[[3]][c(30, 32), 3] <- NA
blocks[[3]][40, ] <- NA
blocks[[4]][42] <- NA
ind_NA <- c(2, 4, 8, 12, 17, 23, 30, 32, 40, 42)

test_that("intersection_list selects the common rows without missing values", {
  blocks_inter <- intersection_list(blocks)
  for (j in seq_along(blocks)) {
    expect_equal(blocks_inter[[j]], subset_rows(blocks[[j]], -ind_NA))
  }
})
test_that("intersection_list raises an error if there is less than 3 subjects
          left after removing missing values", {
  bad_blocks <- blocks
  bad_blocks[[2]][-c(1, 3), 1] <- NA
  expect_error(intersection_list(bad_blocks),
    paste0(
      "Less than 3 subjects have no missing values, ",
      "choose another missing value handling method ",
      "or work on your dataset."
    ),
    fixed = TRUE
  )
})
