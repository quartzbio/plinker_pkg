context('plink lm')

mk_pheno <- function(bo) {
  set.seed(123)
  rnorm(bed_nb_samples(bo))
}

#######################################


.bed_plink_lm <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)
  fam <- bed_fam_df(bo)
  nbind <- bed_nb_samples(bo)

  expect_error(bed_plink_lm(bo, quiet = TRUE), '')

  # must change the phenotype: can not be binary
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
  expect_error(bed_plink_lm(bo, pheno, covars = covars, quiet = TRUE),
    'Error, constant covariates')

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

  # categorical covar
  covars <- fam[, 1:2]
  covars$CATEG1 <- c(rep('toto', nbind - 50), c(rep('titi', 50)))
  set.seed(123)
  covars$CATEG2 <- as.factor(sample(1:10, nrow(covars), replace = TRUE))

  covars <- bed_make_covars(bo, covars)

  res1 <- bed_plink_lm(bo, pheno, covars = covars, quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  ### missing data
  # in phenotype

  pheno3 <- pheno2 <- pheno
  pheno2[11:20] <- NA
  pheno3[11:20] <- -9L # plink default for coding missing
  res1 <- bed_plink_lm(bo, pheno, quiet = TRUE)
  expect_equal(res1$NMISS[1], 89) # no missing

  res2 <- bed_plink_lm(bo, pheno2, quiet = TRUE)
  expect_equal(res2$NMISS[1], 79) # 10 missing

  res3 <- bed_plink_lm(bo, pheno3, quiet = TRUE)
  expect_equal(res3$NMISS[1], 89) # no missing. users do not use -9 but NA

  # in covars
  # add -9L to covars
  covars <- fam[, 1:2]
  covars$COVAR1 <- c(rep(1, nbind - 50), c(rep(2, 50)))
  covars$COVAR1[21:30] <- -9L
  covars2 <- covars
  covars2$COVAR1[21:30] <- NA

  res1 <- bed_plink_lm(bo, pheno, covars = covars, quiet = TRUE)
  expect_equal(res1$NMISS[1], 89) # no missing

  res2 <- bed_plink_lm(bo, pheno, covars = covars2, quiet = TRUE)
  expect_equal(res2$NMISS[1], 79) # 10 missing

  # in both
  res2 <- bed_plink_lm(bo, pheno2, covars = covars2, quiet = TRUE)
  expect_equal(res2$NMISS[1], 69) # 20 missing
}
test_that('bed_plink_lm', .bed_plink_lm())

