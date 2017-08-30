context('missing')

.compute_missing_value_candidate_for_numeric <- function() {
  missing <- plinker:::compute_missing_value_candidate_for_numeric

  expect_error(missing(NULL), 'must be numeric')
  expect_error(missing("toto"), 'must be numeric')

  expect_identical(missing(1), -9L)
  expect_identical(missing(0), -9L)
  expect_identical(missing(0:100), -9L)
  expect_identical(missing(1, candidate = -1), -1L)
  expect_identical(missing(-(0:100)), -101L)
  expect_identical(missing(-257), -9L)
}
test_that('compute_missing_value_candidate_for_numeric',
  .compute_missing_value_candidate_for_numeric())

