context('plink lm')

mk_pheno <- function(bo) {
  set.seed(123)
  rnorm(bed_nb_samples(bo))
}

#######################################


.bed_plink_lm_logistic <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)
  fam <- bed_fam_df(bo)
  nbind <- bed_nb_samples(bo)

  res1 <- bed_plink_lm(bo, logistic = TRUE, quiet = TRUE)
  res2 <- bed_R_lm(bo, logistic = TRUE)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  ### covariates
  fam <- bed_fam_df(bo)
  covars <- fam[, 1:2]

  # constant covar
  covars$CONST <- 1
  expect_error(bed_plink_lm(bo, logistic = TRUE, covars = covars, quiet = TRUE),
    'Error, constant covariates')

  # non constant covars
  covars <- fam[, 1:2]
  covars$COVAR1 <- c(rep(1, nbind - 50), c(rep(2, 50)))
  res1 <- bed_plink_lm(bo, logistic = TRUE, covars = covars, quiet = TRUE)
  res2 <- bed_R_lm(bo, logistic = TRUE, covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  # adding a second covar
  covars$COVAR2 <- fam$SEX
  res1 <- bed_plink_lm(bo, logistic = TRUE, covars = covars, quiet = TRUE)
  res2 <- bed_R_lm(bo, logistic = TRUE, covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  # categorical covar
  covars <- fam[, 1:2]
  covars$CATEG1 <- c(rep('toto', nbind - 50), c(rep('titi', 50)))
  set.seed(123)
  covars$CATEG2 <- as.factor(sample(1:10, nrow(covars), replace = TRUE))

  covars <- bed_make_covars(bo, covars)

  res1 <- bed_plink_lm(bo, logistic = TRUE, covars = covars, quiet = TRUE)
  res2 <- bed_R_lm(bo, logistic = TRUE, covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)
}
test_that('bed_plink_lm_logistic', .bed_plink_lm_logistic())



.bed_plink_lm <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)
  fam <- bed_fam_df(bo)
  nbind <- bed_nb_samples(bo)

  expect_error(bed_plink_lm(bo, quiet = TRUE), 'ERROR running plink')

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
}
test_that('bed_plink_lm', .bed_plink_lm())

