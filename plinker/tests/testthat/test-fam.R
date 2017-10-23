context('fam')

.read_fam <- function() {
  ds <- plinker:::fetch_sample_bed()

  df <- read_fam(plinker:::unprefix_bed(ds)['fam'])

  expect_equal(dim(df), c(89, 6))
  expect_identical(names(df),
    c('FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'SEX', 'PHENO'))

  expect_equal(anyDuplicated(df$IID), 0)

  expect_true(all(unique(df$SEX) %in% 0:2))
  expect_true(all(unique(df$PHENO) %in% 1:2))
}
test_that('read_fam', .read_fam())



.save_fam <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  setup_temp_dir()
  dir.create('dir')
  path <- 'dir/toto.fam'

  fam <- bed_fam_df(bo)

  save_fam(fam, path)

  fam2 <- read_fam(path)
  expect_identical(fam2, fam)
}
test_that('save_fam', .save_fam())



.merge_df_with_fam <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  fam_df <- bed_fam_df(bo)
  nb <- nrow(fam_df)

  ################################### by FID,IID ###

  # with itself
  mdf <- merge_df_with_fam(fam_df, fam_df)
  expect_identical(mdf[, 1:6], fam_df)

  # itself, but reversed
  mdf <- merge_df_with_fam(fam_df, fam_df[nb:1, ])
  expect_identical(mdf[, 1:6], fam_df)

  # subset
  expect_error(merge_df_with_fam(fam_df, fam_df[1:10, ]), 'bad merge')

  # superset
  df <- rbind(fam_df[1,], fam_df)
  df[1, 'IID'] <- 'somethingunique'
  mdf <- merge_df_with_fam(fam_df, df)
  expect_identical(mdf[, 1:6], fam_df)

  # repeated rows
  df <- rbind(fam_df[1,], fam_df)
  expect_error(merge_df_with_fam(fam_df, df), 'bad merge')

  # rm_merge_id
  df <- fam_df
  df$NB <- 1:nrow(df)
  mdf1 <- merge_df_with_fam(fam_df, df, rm_merge_id = TRUE)
  mdf2 <- merge_df_with_fam(fam_df, df, rm_merge_id = FALSE)
  # merge_id not used for merging
  expect_identical(mdf1, mdf2)


  ############ by single ID ##################
  ## ignore_fid = FALSE
  ids <- bed_sample_IDs(bo, ignore_fid = FALSE)

  df <- data.frame(MYID = rev(ids), VALUE = seq_along(ids),
    stringsAsFactors = FALSE)
  expect_error(merge_df_with_fam(fam_df, df, merge_id = 'SUBJID', ignore_fid = FALSE),
    'error, bad "merge_id"')
  expect_error(merge_df_with_fam(fam_df, df, merge_id = 'MYID',ignore_fid = FALSE),
    NA)
  expect_error(merge_df_with_fam(fam_df, df,ignore_fid = FALSE),
    NA)


  df <- data.frame(SUBJID = rev(ids), VALUE = seq_along(ids),
     stringsAsFactors = FALSE)

  mdf <- merge_df_with_fam(fam_df, df, ignore_fid = FALSE)

  expect_identical(mdf[, -7], fam_df)
  expect_identical(mdf$VALUE, rev(seq_along(ids)))

  # rm_merge_id
  mdf2 <- merge_df_with_fam(fam_df, df, ignore_fid = FALSE, rm_merge_id = FALSE)
  expect_identical(mdf2$SUBJID, ids)
  mdf2$SUBJID <- NULL
  expect_identical(mdf2, mdf)

  # subset
  expect_error(merge_df_with_fam(fam_df,  df[1:10, ],
      ignore_fid = FALSE), 'bad merge')

  # superset
  df2 <- rbind(df[1, ], df)
  df2[1, 1] <- 'COCORICO'
  mdf2 <- merge_df_with_fam(fam_df, df2, ignore_fid = FALSE)
  expect_identical(mdf2, mdf)

  # repeated rows
  df2 <- rbind(df[1, ], df)
  expect_error(merge_df_with_fam(fam_df,  df2,
      ignore_fid = FALSE), 'bad merge')

  ## ignore_fid = TRUE
  ids <- bed_sample_IDs(bo, ignore_fid = TRUE)
  df <- data.frame(SUBJID = rev(ids), VALUE = seq_along(ids),
    stringsAsFactors = FALSE)

  expect_error(merge_df_with_fam(fam_df, df, merge_id = 'MYID', ignore_fid = TRUE),
    'error, bad "merge_id"')

  mdf <- merge_df_with_fam(fam_df, df, ignore_fid = TRUE)
  expect_identical(mdf[, -7], fam_df)
  expect_identical(mdf$VALUE, rev(seq_along(ids)))

  # mismatch in ignore_fid
  expect_error(merge_df_with_fam(fam_df,  df, ignore_fid = FALSE), 'bad merge')

}
test_that('merge_df_with_fam', .merge_df_with_fam())
