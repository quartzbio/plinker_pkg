context('bed object')


.bed_open <- function() {
  prefix <- plinker:::fetch_sample_bed()

  bo <- bed_open(prefix)

  expect_is(bo, 'plinker_bed')
  expect_equal(bo$nb_snps, 17)
  expect_equal(nrow(bo$fam_df), 89)
  expect_false(bo$ignore_fid)
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
  expect_null(bed_sample_idx(bo))
  expect_identical(bed_fam_df(bo), read_fam(bo$fam))
  expect_identical(bed_bim_df(bo), read_bim(bo$bim))

  expect_identical(bed_snp_IDs(bo), bed_bim_df(bo)$SNP)

  expect_false(bed_ignore_fid(bo))

  bim <- bed_bim_df(bo)
  expect_identical(bed_allele1(bo), bim$A1)
  expect_identical(bed_allele2(bo), bim$A2)

  expect_identical(bed_allele_higher(bo), pmax(bim$A2, bim$A1))


  bo2 <- bed_subset(bo, snp_idx = 1)
  expect_identical(bed_allele1(bo2), bim$A1[1])
  expect_identical(bed_allele2(bo2), bim$A2[1])
}
test_that('accessors', .accessors())



.print.plinker_bed <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  expect_output(print(bo), 'PLINK BED')
  expect_output(print(bo), 'Not yet loaded')
  expect_output(print(bo), 'sample IDs: FID_IID')

  df <- bed_bim_df(bo) # trigger the bim_df in cache
  expect_output(print(bo), 'Loaded')

  ### sample_annot
  ids <- bed_sample_IDs(bo)
  df <- data.frame(MERGE_ID = rev(ids), SUBJID = as.character(seq_along(ids)),
    stringsAsFactors = FALSE)

  bo2 <- bed_set_sample_annot(bo, df, id = 'SUBJID')
  expect_output(print(bo2), 'sample annotations: custom ID="SUBJID", 1 var(s)',
    fixed = TRUE)

  bo2 <- bed_set_sample_annot(bo, df)
  expect_output(print(bo2), 'sample annotations: 1 var(s)', fixed = TRUE)
}
test_that('print.plinker_bed', .print.plinker_bed())



.as.data.frame.plinker_bed <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  bo2 <- bed_subset(bo, snp_idx = 5:13, sample_idx = 24:67)
  expect_identical(as.data.frame(bo2),
    plinker:::bed_convert_genotypes_to_data_frame(bo2))
}
test_that('as.data.frame.plinker_bed', .as.data.frame.plinker_bed())



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



.bed_sample_IDs <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  fam_df <- bed_fam_df(bo)

  ids <- bed_sample_IDs(bo)
  expect_identical(ids, paste0(fam_df$FID, '_', fam_df$IID))

  ids <- bed_sample_IDs(bo, ignore_fid = TRUE)
  expect_identical(ids, fam_df$IID)

  bo2 <- bed_subset_samples_by_idx(bo, 5:20)
  ids <- bed_sample_IDs(bo2, ignore_fid = TRUE)
  expect_identical(ids, fam_df$IID[5:20])

  ids <- bed_sample_IDs(bo2, ignore_fid = FALSE)
  expect_identical(ids, paste0(fam_df$FID, '_', fam_df$IID)[5:20])
}
test_that('bed_sample_IDs', .bed_sample_IDs())

