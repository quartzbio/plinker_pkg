context('phenotype')


.bed_phenotype_from_df <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  fam <- bed_fam_df(bo)
  nb <- nrow(fam)

  pheno <- bed_phenotype_from_df(bo, fam, 'PHENO')
  expect_identical(pheno, fam$PHENO)

  ids <- bed_sample_IDs(bo)
  df <- data.frame(SUBJID = rev(ids), VALUE = nchar(ids),
    stringsAsFactors = FALSE)

  df[nb:10, 'VALUE'] <- NA
  pheno <- bed_phenotype_from_df(bo, df, 'VALUE', id_var = 'SUBJID')

  expect_error(make_phenotype_from_vector(pheno, missing_phenotype = 0L),
    'must be negative integer')

  pheno <- make_phenotype_from_vector(pheno, missing_phenotype = -1L)

  expect_true(all(tail(pheno, 9) != -1))
  expect_true(all(head(pheno, nb - 10) == -1))
}
test_that('bed_phenotype_from_df', .bed_phenotype_from_df())



.make_phenotype_from_vector <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  fam <- bed_fam_df(bo)

  pheno <- make_phenotype_from_vector(fam$PHENO)
  expect_identical(pheno, fam$PHENO)
  expect_true(is.integer(pheno))

  v <- fam$PHENO
  v[1] <- v[17] <- NA
  pheno <- make_phenotype_from_vector(v)
  expect_equal(which(pheno == -9), c(1, 17))

  # safety net
  v2 <- v
  v2[23] <- -9L

  expect_error(make_phenotype_from_vector(v2), 'missing_phenotype')
  # but no error if no NA
  v2[1] <- v2[17] <- -9L
  expect_identical(make_phenotype_from_vector(v2), v2)

  pheno <- make_phenotype_from_vector(v, missing_phenotype = -2L)
  expect_equal(which(pheno == -2), c(1, 17))

  expect_error(make_phenotype_from_vector(v, missing_phenotype = 0),
    'bad missing_phenotype')
  expect_error(make_phenotype_from_vector(v, missing_phenotype = NA),
    'bad missing_phenotype')
  expect_error(make_phenotype_from_vector(v, missing_phenotype = 0),
    'bad missing_phenotype')
  expect_error(make_phenotype_from_vector(v, missing_phenotype = "0"),
    'bad missing_phenotype')
}
test_that('make_phenotype_from_vector', .make_phenotype_from_vector())

