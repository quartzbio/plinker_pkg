context('bed object')

.new_bed <- function() {
  new_bed <- plinker:::new_bed

  fns <- as.list(plinker:::unprefix_bed(plinker:::fetch_sample_bed()))

  ### edge cases
  expect_error(new_bed(bed = fns$bed, bim = fns$bim, fam = 'foo'), 'fam file')
  expect_error(new_bed(bed = 'foo', bim = fns$bim, fam = fns$fam), 'bed file')
  expect_error(new_bed(bed = fns$bed, bim = 'bar', fam = fns$fam), 'bim file')

  bed <- new_bed(bed = fns$bed, bim = fns$bim, fam = fns$fam)
  expect_is(bed, 'plinker_bed')

  expect_identical(bed[names(fns)], fns)
}
test_that('new_bed', .new_bed())

