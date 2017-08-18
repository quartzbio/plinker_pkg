context('phenotype')

.bed_phenotype_from_vector <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  fam <- bed_fam_df(bo)

  pheno <- bed_phenotype_from_vector(bo, fam$PHENO)
  expect_identical(pheno, fam$PHENO)
  expect_true(is.integer(pheno))

  v <- fam$PHENO
  v[1] <- v[17] <- NA
  pheno <- bed_phenotype_from_vector(bo, v)
  expect_equal(which(pheno == -9), c(1, 17))

  # safety net
  v[23] <- -9L
  expect_error(bed_phenotype_from_vector(bo, v), 'missing_phenotype')

  pheno <- bed_phenotype_from_vector(bo, v, missing_phenotype = 0L)
  expect_equal(which(pheno == 0), c(1, 17))

  expect_error(bed_phenotype_from_vector(bo, v, missing_phenotype = 0),
    'bad missing_phenotype')
  expect_error(bed_phenotype_from_vector(bo, v, missing_phenotype = NA),
    'bad missing_phenotype')
  expect_error(bed_phenotype_from_vector(bo, v, missing_phenotype = 0),
    'bad missing_phenotype')
  expect_error(bed_phenotype_from_vector(bo, v, missing_phenotype = "0"),
    'bad missing_phenotype')
}
test_that('bed_phenotype_from_vector', .bed_phenotype_from_vector())

