context('convert')

.bed_convert_genotypes_to_data_frame <- function() {
  recode_genotypes <- plinker:::recode_genotypes

  bo <- bed_open(plinker:::fetch_sample_bed())

  bo1 <- bed_subset(bo, snp_idx = 2, sample_idx = 1:20)

  df <- bed_convert_genotypes_to_data_frame(bo1)
  expect_is(df, 'data.frame')
  expect_equal(nrow(df), bed_nb_samples(bo1) * bed_nb_snps(bo1))

  expect_equivalent(unique(df[, 1:6]), bed_bim_df(bo1))
  genos <- unname(bed_genotypes(bo1)[, 1])
  expect_identical(df$GENO_INT, genos)
  expect_identical(df$GENO_STR, unname(bed_genotypes_as_strings(bo1)[, 1]))

  expect_identical(df$DOMINANT, recode_genotypes(genos, 'dominant'))
  expect_identical(df$RECESSIVE, recode_genotypes(genos, 'recessive'))

  ### scalar df
  bo1 <- bed_subset(bo, snp_idx = 10, sample_idx = 80)
  df <- bed_convert_genotypes_to_data_frame(bo1)
  expect_equal(nrow(df), 1)

  genos <- unname(bed_genotypes(bo1)[, 1])
  expect_identical(df$GENO_INT, genos)
  expect_identical(df$GENO_STR, unname(bed_genotypes_as_strings(bo1)[, 1]))

  ### all dataset
  df <- bed_convert_genotypes_to_data_frame(bo)
  expect_equal(nrow(df), bed_nb_samples(bo) * bed_nb_snps(bo))

  expect_equivalent(unique(df[, 1:6]), bed_bim_df(bo))
  genos <- as.integer(bed_genotypes(bo))
  expect_identical(df$GENO_INT, genos)
  expect_identical(df$GENO_STR, as.character(bed_genotypes_as_strings(bo)))
}
test_that('bed_convert_genotypes_to_data_frame',
  .bed_convert_genotypes_to_data_frame())

