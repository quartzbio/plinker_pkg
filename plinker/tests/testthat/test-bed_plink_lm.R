context('plink lm')




.bed_plink_lm <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)
  fam <- bed_fam_df(bo)
  nbind <- bed_nb_samples(bo)

  expect_error(bed_plink_lm(bo), '')

  # must change the phenotype: can not be binary
  mk_pheno <- function(bo) {
      set.seed(123)
      rnorm(bed_nb_samples(bo))
  }
  pheno <- mk_pheno(bo)
  res1 <- bed_plink_lm(bo, pheno, quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  ### covariates
  fam <- bed_fam_df(bo)
  covars <- fam[, 1:2]

  # constant covar
  covars$CONST <- 1
  bed_plink_lm(bo, pheno)
  res1 <- bed_plink_lm(bo, pheno, covars = covars, quiet = TRUE)
  # ==> plink output all NA
  expect_true(all(is.na(res1[, 7:9])))

  # non constant covars
  covars <- fam[, 1:2]
  covars$COVAR1 <- c(rep(1, nbind - 50), c(rep(2, 50)))
  res1 <- bed_plink_lm(bo, pheno, covars = covars, quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  # adding a second covar
  covars$COVAR2 <- fam$SEX
  res1 <- bed_plink_lm(bo, pheno, covars = covars, quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)
  browser()
   # categorical covar
   covars <- fam[, 1:2]
   covars$CATEG <- c(rep('toto', nbind - 50), c(rep('titi', 50)))

   res1 <- bed_plink_lm(bo, pheno, covars = bed_make_covars(bo, covars),
     quiet = TRUE)

   res2 <- bed_R_lm(bo, covars = covars)
   covars$CATEG <- as.factor(covars$CATEG)
   res3 <- bed_R_lm(bo, covars = covars)


}
test_that('bed_plink_lm', .bed_plink_lm())

