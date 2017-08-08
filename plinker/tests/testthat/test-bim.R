context('bim')

.read_bim <- function() {
  ds <- plinker:::fetch_sample_bed()

  df <- read_bim(plinker:::unprefix_bed(ds)['bim'])

  expect_equal(dim(df), c(17, 6))
  expect_identical(names(df),
    c('CHR', 'SNPID', 'MORGANS', 'POS', 'A1', 'A2'))

  expect_identical(unique(df$CHR), "8")
  expect_equal(anyDuplicated(df$SNPID), 0)

  expect_equal(anyDuplicated(df$POS), 0)

  expect_true(all(unique(df$A1) %in% c('A', 'C', 'G', 'T')))
  expect_true(all(unique(df$A2) %in% c('A', 'C', 'G', 'T')))
}
test_that('read_bim', .read_bim())

