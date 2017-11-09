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

  ### allele order: inversed OR, opposite stat, same pvalue
  ### PB: bug in PLINK, lots of NA
  res1 <- bed_plink_lm(bo, logistic = TRUE, quiet = TRUE)
  res2 <- bed_plink_lm(bo, logistic = TRUE, quiet = TRUE,
    lexicographic_allele_order = TRUE)
  res3 <- bed_R_lm(bo, logistic = TRUE, lexicographic_allele_order = TRUE)

  non_missing <- which(!is.na(res2$P))
  expect_equal(res3[non_missing, ], res2[non_missing, ], tolerance = 2e-4)
  # non swapped allele
  expect_identical(res2[3, ], res1[3, ])

  # swapped allele
  expect_equal(res2[14, 'OR'], 1/res1[14, 'OR'], tolerance = 1e-4)
  expect_equal(res2[14, 'STAT'], -res1[14, 'STAT'], tolerance = 1e-4)
  expect_equal(res2[14, 'P'], res1[14, 'P'], tolerance = 1e-4)
}
test_that('bed_plink_lm_logistic', .bed_plink_lm_logistic())



.bed_plink_lm_recessive <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)
  fam <- bed_fam_df(bo)
  nbind <- bed_nb_samples(bo)

  # must change the phenotype: can not be binary
  pheno <- mk_pheno(bo)
  res1 <- bed_plink_lm(bo, pheno, model = 'recessive', quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, model = 'recessive')

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  ### covariates
  fam <- bed_fam_df(bo)
  covars <- fam[, 1:2]

  # non constant covars
  covars <- fam[, 1:2]
  covars$COVAR1 <- c(rep(1, nbind - 50), c(rep(2, 50)))
  res1 <- bed_plink_lm(bo, pheno, covars = covars, model = 'recessive', quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, model = 'recessive', covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])

  # ==> more NA in PLINK
  na1 <- which(is.na(res1$P))
  na2 <- which(is.na(res2$P))
  expect_length(na1, 8)
  expect_length(na2, 4)
  expect_true(all(na2 %in% na1))

  # equal on non NA
  expect_equal(res2[-na1, 7:9], res1[-na1, 7:9], tolerance = 1e-3)

  # adding a second covar
  covars$COVAR2 <- fam$SEX
  res1 <- bed_plink_lm(bo, pheno, model = 'recessive', covars = covars, quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, model = 'recessive', covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  na1 <- which(is.na(res1$P))
  expect_equal(res2[-na1, 7:9], res1[-na1, 7:9], tolerance = 1e-3)

  # categorical covar
  covars <- fam[, 1:2]
  covars$CATEG1 <- c(rep('toto', nbind - 50), c(rep('titi', 50)))
  set.seed(123)
  covars$CATEG2 <- as.factor(sample(1:10, nrow(covars), replace = TRUE))

  covars <- bed_make_covars(bo, covars)

  res1 <- bed_plink_lm(bo, pheno, covars = covars, model = 'recessive', quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, covars = covars, model = 'recessive')

  expect_identical(res2[, 1:6], res1[, 1:6])
  na1 <- which(is.na(res1$P))
  expect_equal(res2[-na1, 7:9], res1[-na1, 7:9], tolerance = 1e-3)
}
test_that('bed_plink_lm_recessive', .bed_plink_lm_recessive())



.bed_plink_lm_dominant <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)
  fam <- bed_fam_df(bo)
  nbind <- bed_nb_samples(bo)

  # must change the phenotype: can not be binary
  pheno <- mk_pheno(bo)
  res1 <- bed_plink_lm(bo, pheno, model = 'dominant', quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, model = 'dominant')

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  ### covariates
  fam <- bed_fam_df(bo)
  covars <- fam[, 1:2]

  # non constant covars
  covars <- fam[, 1:2]
  covars$COVAR1 <- c(rep(1, nbind - 50), c(rep(2, 50)))
  res1 <- bed_plink_lm(bo, pheno, covars = covars, model = 'dominant', quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, model = 'dominant', covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

  # adding a second covar
  covars$COVAR2 <- fam$SEX
  res1 <- bed_plink_lm(bo, pheno, model = 'dominant', covars = covars, quiet = TRUE)
  res2 <- bed_R_lm(bo, pheno, model = 'dominant', covars = covars)

  expect_identical(res2[, 1:6], res1[, 1:6])
  expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)

   # categorical covar
   covars <- fam[, 1:2]
   covars$CATEG1 <- c(rep('toto', nbind - 50), c(rep('titi', 50)))
   set.seed(123)
   covars$CATEG2 <- as.factor(sample(1:10, nrow(covars), replace = TRUE))

   covars <- bed_make_covars(bo, covars)

   res1 <- bed_plink_lm(bo, pheno, covars = covars, model = 'dominant', quiet = TRUE)
   res2 <- bed_R_lm(bo, pheno, covars = covars, model = 'dominant')

   expect_identical(res2[, 1:6], res1[, 1:6])
   expect_equal(res2[, 7:9], res1[, 7:9], tolerance = 1e-3)


}
test_that('bed_plink_lm_dominant', .bed_plink_lm_dominant())



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

  ### ordering
  bo2 <- bed_subset(bo, snp_IDs = rev(bed_snp_IDs(bo)))
  res <- bed_plink_lm(bo2, pheno, covars = covars, quiet = TRUE)
  ref <- bed_plink_lm(bo, pheno, covars = covars, quiet = TRUE)

  expect_false(identical(res, ref))
  res2 <- plinker:::reorder_plink_snp_output(res, bed_snp_IDs(bo))
  expect_identical(res2, ref)

  res3 <- bed_R_lm(bo2, pheno, covars = covars)
  expect_identical(res3[, 1:6], res[, 1:6])


  ### allele order: negate beta and stat, same pvalue
  res1 <- bed_plink_lm(bo, pheno, quiet = TRUE)
  res2 <- bed_plink_lm(bo, pheno, lexicographic_allele_order = TRUE, quiet = TRUE)
  res3 <- bed_R_lm(bo, pheno, lexicographic_allele_order = TRUE)
  expect_equal(res3, res2, tolerance = 2e-4)

  expect_identical(res1[1, ], res2[1, ])
  expect_equal(res2$P[2], res1$P[2])
  expect_equal(res2$BETA[2], -res1$BETA[2])
  expect_equal(res2$STAT[2], -res1$STAT[2])


  ### snp annot
  bim <- bed_bim_df(bo)
  annot <- data.frame(
    SNP = bim$SNP,
    ID = paste0('ID_', bim$SNP),
    TOTO = 1,
    stringsAsFactors = FALSE)

  bo2 <- bed_set_snp_annot(bo, annot, 'ID')
  res <- bed_plink_lm(bo, pheno, quiet = TRUE)
  res2 <- bed_plink_lm(bo2, pheno, quiet = TRUE)

  expect_identical(res2[, -3], res)
  expect_identical(paste0('ID_', res2$SNP), res2$ID)
}
test_that('bed_plink_lm', .bed_plink_lm())

