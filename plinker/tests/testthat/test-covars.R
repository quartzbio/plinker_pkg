context('covariates')


.dummify_df_vars <- function() {
  dummify_df_vars <- plinker:::dummify_df_vars

#  expect_error(dummify_factor(1), 'give a factor')

  mat <- dummify_df_vars(iris[, 'Species', drop = FALSE])
  expect_is(mat, 'matrix')
  expect_identical(colnames(mat), paste0('Species', levels(iris$Species)[-1]))
  expect_equal(nrow(mat), nrow(iris))

  mat <- dummify_df_vars(iris)
  expect_identical(colnames(mat),
    c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width",
      paste0('Species', levels(iris$Species)[-1])))

  ###
  df <- data.frame(
    INT = 1:3,
    NUM = as.numeric(11:13),
    STR = LETTERS[1:3],
    FACTOR = as.factor(LETTERS[1:3]),
    STR2 = letters[1:3],
    stringsAsFactors = FALSE
  )
  mat <- dummify_df_vars(df)
  expect_identical(colnames(mat),
    c("INT", "NUM", "STRB", "STRC", "FACTORB", "FACTORC", "STR2b", "STR2c"))
}
test_that('dummify_df_vars', .dummify_df_vars())



.check_plink_covars <- function() {
  check_plink_covars <- plinker:::check_plink_covars

  bo <- bed_open(plinker:::fetch_sample_bed())
  fam <- bed_fam_df(bo)

  expect_error(check_plink_covars(fam, fam),
    'Error, categorical covariates: FATHER_ID, MOTHER_ID')

  expect_error(check_plink_covars(iris, fam), 'missing required')

  expect_error(check_plink_covars(head(fam), fam), 'do not match')

  ### POSITIVE CONTROL
  df <- fam[, c('FID', 'IID', 'SEX', 'PHENO')]
  expect_error(check_plink_covars(df, fam), NA)

  # does not rely on rownames
  rownames(df) <- rev(rownames(fam))
  expect_error(check_plink_covars(df, fam), NA)

  ### constant covars
  covars <- fam[, 1:2]
  covars$COVAR1 <- 1:nrow(covars)
  covars$CONST1 <- 1
  covars$COVAR2 <- rev(covars$COVAR1)
  covars$CONST2 <- 2
  covars$CONST2[1] <- NA

  expect_error(check_plink_covars(covars, fam),
    'constant covariates: CONST1, CONST2')

  ### categorical variables
  covars <- fam[, 1:2]
  covars$STR <- 'toto'
  covars$STR[2] <- 'titi'
  covars$COVAR1 <- 1:nrow(covars)
  covars$FACTOR <- as.factor(1:nrow(covars))
  expect_error(check_plink_covars(covars, fam),
    'categorical covariates: STR, FACTOR')
}
test_that('check_plink_covars', .check_plink_covars())




.bed_make_covars <- function() {
  check_plink_covars <- plinker:::check_plink_covars

  bo <- bed_open(plinker:::fetch_sample_bed())

  fam <- bed_fam_df(bo)

  ids <- bed_sample_IDs(bo)


  ##### subjid

  covars <- data.frame(
    SUBJID = rev(ids),
    TOTO = nchar(rev(ids)),
    stringsAsFactors = FALSE)

  df <- bed_make_covars(bo, covars)

  expect_identical(df[, 1:2], fam[, 1:2])
  expect_equal(ncol(df), 3)
  expect_equal(df$TOTO, nchar(df$FID) + nchar(df$IID) + 1)

  ##### test 2 - id var
  covars$MYID <- covars$SUBJID
  covars$SUBJID <- NULL
  expect_error(bed_make_covars(bo, covars, 'SUBJID'), 'bad "merge_id"')

  df2 <- bed_make_covars(bo, covars, merge_id = 'MYID')
  expect_identical(df2, df)

  ##### test 3 - FID, IID
  covars <- fam[, 1:2]
  covars$COVAR1 <- 1:nrow(covars)
  covars$COVAR2 <- rev(covars$COVAR1)
  covars <- covars[nrow(covars):1, ]

  df <- bed_make_covars(bo, covars)
  expect_identical(df[, 1:2], fam[, 1:2])
  expect_identical(names(df), c("FID", "IID", "COVAR1", "COVAR2"))
  expect_equal(df$COVAR1, 1:nrow(covars))
  expect_equal(df$COVAR2, rev(df$COVAR1))

  ### with categorical var
  covars <- fam[, 1:2]
  covars$CATEG <- substr(covars$FID, 1, 3)

  df <- bed_make_covars(bo, covars)
  expect_identical(df[, 1:2], fam[, 1:2])
  fct <- as.factor(covars$CATEG)
  expect_identical(names(df)[-(1:2)], paste0('CATEG', levels(fct)[-1]))




}
test_that('bed_make_covars', .bed_make_covars())

