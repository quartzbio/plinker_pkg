context('utils')

.split_sorted_ints_by_blocks <- function() {
  split <- plinker:::split_sorted_ints_by_blocks

  expect_null(split(NULL))

  expect_identical(split(1L), cbind(1L, 1L))

  expect_identical(split(1:2), cbind(1L, 2L))

  expect_identical(split(c(3L, 10:20, 25L, 30:100)),
    cbind(
      c(3L, 10L, 25L, 30L),
      c(1L, 11L, 1L, 71L)
      ))
}
test_that('split_sorted_ints_by_blocks', .split_sorted_ints_by_blocks())

