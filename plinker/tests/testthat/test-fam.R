context('fam')

.read_fam <- function() {
  ds <- plinker:::fetch_sample_bed()

  df <- read_fam(plinker:::unprefix_bed(ds)['fam'])

  expect_equal(dim(df), c(89, 6))
  expect_identical(names(df),
    c('FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'SEX', 'PHENO'))
  expect_equal(anyDuplicated(df$FID), 0)
  expect_equal(anyDuplicated(df$IID), 0)

  expect_true(all(unique(df$SEX) %in% 1:2))
  expect_true(all(unique(df$PHENO) %in% 1:2))
}
test_that('read_fam', .read_fam())

