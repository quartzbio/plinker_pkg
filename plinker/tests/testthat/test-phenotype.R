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
  expect_identical(pheno, rev(df$VALUE))
}
test_that('bed_phenotype_from_df', .bed_phenotype_from_df())


