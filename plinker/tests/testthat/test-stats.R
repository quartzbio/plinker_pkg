context('R stats')

.bed_R_lm <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  bo1 <- bed_subset(bo, snp_idx = 10:12)
  nbind <- bed_nb_samples(bo1)

  res <- bed_R_lm(bo1)
  expect_is(res, 'data.frame')
  expect_equal(nrow(res), bed_nb_snps(bo1))
  expect_identical(names(res),
    c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P"))

  # reproducible
  expect_identical(bed_R_lm(bo1), res)

  ### phenotype
  pheno <- 1:nbind
  bo1 <- bed_subset(bo, snp_idx = 10:12)
  res2 <- bed_R_lm(bo1, pheno)
  expect_false(identical(res2, res))

  ### covars
  # quantitative
  fam <- bed_fam_df(bo1)
  covars <- fam[, 1:2]

  # constant covar
  covars$CONST <- 1
  res2 <- bed_R_lm(bo1, covars = covars)

  expect_equivalent(res2[res2$TEST == 'ADD', ], res)

  # non constant covars
  covars <- fam[, 1:2]
  covars$COVAR1 <- c(rep(1, nbind - 50), c(rep(2, 50)))
  res2 <- bed_R_lm(bo1, covars = covars)
  expect_true(abs(res2$P[1] - res$P[1]) > 0.05)

  # categorical covar
  covars <- fam[, 1:2]
  covars$CATEG <- c(rep('toto', nbind - 50), c(rep('titi', 50)))

  res2 <- bed_R_lm(bo1, covars = covars)
  covars$CATEG <- as.factor(covars$CATEG)
  res3 <- bed_R_lm(bo1, covars = covars)
  expect_identical(res3, res2)

}
test_that('bed_R_lm', .bed_R_lm())

