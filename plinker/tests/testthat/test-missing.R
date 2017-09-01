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



.compute_missing_value_candidate <- function() {
  missing <- plinker:::compute_missing_value_candidate

  bo <- bed_open(plinker:::fetch_sample_bed())
  fam <- bed_fam_df(bo)

  ### constant covars
  covars <- data.frame(COVAR1 = 1:nrow(fam))
  covars$CONST1 <- 1
  covars$COVAR2 <- rev(covars$COVAR1)
  covars$CONST2 <- 2

  # no missing yet
  expect_identical(missing(fam$PHENO, covars), -9L)
  expect_identical(missing(fam$PHENO), -9L)
  expect_identical(missing(fam$PHENO, covars, candidate = -3L), -3L)
  expect_identical(missing(fam$PHENO, candidate = -3L), -3L)

  # add -9L to covars
  covars$COVAR2[3] <- -9
  expect_identical(missing(fam$PHENO, covars), -10L)
}
test_that('compute_missing_value_candidate',
  .compute_missing_value_candidate())
