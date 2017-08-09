context('bim')


#.read_bim_snpids_by_block <- function() {
#  read <- plinker:::read_bim_snpids_by_block
#  bim <- plinker:::unprefix_bed(plinker:::fetch_sample_bed())['bim']
#  split_idx <- plinker:::split_sorted_ints_by_blocks
#
#  df0 <- read_bim(bim)
#
#  blks <- split_idx(c(1:5, 10:nrow(df0)))
#  df <- read(bim, blks)
#
#
#
#}
#test_that('read_bim_snpids_by_block', .read_bim_snpids_by_block())



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

