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
  pheno <- bed_phenotype_from_df(bo, df, 'VALUE', id_var = 'SUBJID',
    missing_phenotype = 0L)

  expect_true(all(tail(pheno, 9) != 0))
  expect_true(all(head(pheno, nb - 10) == 0))
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
  v[23] <- -9L
  expect_error(make_phenotype_from_vector(v), 'missing_phenotype')

  pheno <- make_phenotype_from_vector(v, missing_phenotype = 0L)
  expect_equal(which(pheno == 0), c(1, 17))

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

