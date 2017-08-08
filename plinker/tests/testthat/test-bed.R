context('bed object')


.bed_open <- function() {
  prefix <- plinker:::fetch_sample_bed()

  obj <- bed_open(prefix)

  expect_is(obj, 'plinker_bed')
  expect_equal(obj$nb_snps, 17)
  expect_equal(nrow(obj$fam_df), 89)
}
test_that('bed_open', .bed_open())



.new_bed <- function() {
  new_bed <- plinker:::new_bed

  fns <- as.list(plinker:::unprefix_bed(plinker:::fetch_sample_bed()))

  ### edge cases
  expect_error(new_bed(bed = fns$bed, bim = fns$bim, fam = 'foo'), 'fam file')
  expect_error(new_bed(bed = 'foo', bim = fns$bim, fam = fns$fam), 'bed file')
  expect_error(new_bed(bed = fns$bed, bim = 'bar', fam = fns$fam), 'bim file')

  bed <- new_bed(bed = fns$bed, bim = fns$bim, fam = fns$fam,
    fam_df = NULL, nb_snps = 2)
  expect_is(bed, 'plinker_bed')

  expect_identical(bed[names(fns)], fns)
}
test_that('new_bed', .new_bed())



.accessors <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  expect_equal(bed_nb_samples(bo), 89)
  expect_equal(bed_nb_snps(bo), 17)
  expect_null(bed_snp_idx(bo))
  expect_identical(bed_fam_df(bo), read_fam(bo$fam))
  expect_identical(bed_bim_df(bo), read_bim(bo$bim))
}
test_that('accessors', .accessors())



.print.plinker_bed <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  expect_output(print(bo), 'PLINK BED')
  expect_output(print(bo), 'Not yet loaded')
}
test_that('print.plinker_bed', .print.plinker_bed())



.infer_nb_snps <- function() {
  infer_nb_snps <- plinker:::infer_nb_snps

  fns <- as.list(plinker:::unprefix_bed(plinker:::fetch_sample_bed()))
  fam_df <- read_fam(fns$fam)
  nb_samples <- nrow(fam_df)

  nb_snps <- infer_nb_snps(fns$bed, nb_samples)

  bim_df <- read_bim(fns$bim)
  expect_equal(nb_snps, nrow(bim_df))
}
test_that('infer_nb_snps', .infer_nb_snps())
