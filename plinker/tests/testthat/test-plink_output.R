context('plink output')

.annotate_plink_snp_output <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  annotate_plink_snp_output <- plinker:::annotate_plink_snp_output

  bim <- bed_bim_df(bo)
  ids <- bed_snp_IDs(bo)

  # create a PLINK-like output
  # duplicates, missing and some in reverse order
  df <- bim[1:5, ]
  df <- rbind(df, df[nrow(df):1, ])

  df$A1 <- df$A2 <- df$MORGANS <- NULL
  df$RESULT <- 1:nrow(df)

  annot <- data.frame(
    SNP = bim$SNP,
    ID = paste0('ID_', bim$SNP),
    TOTO = 1,
    stringsAsFactors = FALSE)

  expect_error(annotate_plink_snp_output(df, annot, 'ARF'), 'bad param annot_id')
  expect_error(annotate_plink_snp_output(df, annot[, -1], 'ID'), 'mandatory cols')

  df2 <- annotate_plink_snp_output(df, annot, 'ID')

  expect_identical(df2[, -3], df)
  # check ID mapping
  expect_identical(paste0('ID_', df2$SNP), df2$ID)
}
test_that('annotate_plink_snp_output', .annotate_plink_snp_output())



.reorder_plink_sample_output <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  ignore_fid = bed_ignore_fid(bo)
  reorder_plink_sample_output <- plinker:::reorder_plink_sample_output

  fam <- bed_fam_df(bo)
  ids <- bed_sample_IDs(bo)

  ### edge cases
  expect_error(reorder_plink_sample_output(iris, ids), 'missing FID')
  df <- fam
  df$IID <- paste0('toto', df$IID)
  expect_error(reorder_plink_sample_output(df, ids), 'bad sample ids')

  expect_identical(reorder_plink_sample_output(fam, NULL, ignore_fid), fam)

  df2 <- reorder_plink_sample_output(fam, ids, ignore_fid)
  expect_identical(df2, fam)

  ### ignore_fid
  expect_error(reorder_plink_sample_output(fam, ids, !ignore_fid),
     'bad sample ids')

  df <- fam[nrow(fam):1, ]
  df2 <- reorder_plink_sample_output(df, ids, ignore_fid)
  expect_identical(df2, fam)

  df3 <- reorder_plink_sample_output(fam, rev(ids), ignore_fid)
  expect_equivalent(df3, df)

  ### smaller sets
  df2 <- reorder_plink_sample_output(fam[1:10, ], rev(ids), ignore_fid)
  expect_equivalent(df2, fam[10:1, ])

  ### output may have several values by SNP (long format)
  df <- fam
  idx <- rep(1:nrow(fam), each = 4)
  df <- fam[idx, ]
  df$IDX <- 1:nrow(df)
  rownames(df) <- NULL

  df2 <- reorder_plink_sample_output(df, ids, ignore_fid)
  expect_identical(df2, df)

  df2 <- reorder_plink_sample_output(df, rev(ids), ignore_fid)
  ids2 <- compute_sample_IDs(df2, ignore_fid)
  expect_identical(unique(ids2), rev(ids))
}
test_that('reorder_plink_sample_output', .reorder_plink_sample_output())



.annotate_plink_sample_output <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  ignore_fid = bed_ignore_fid(bo)
  annotate_plink_sample_output <- plinker:::annotate_plink_sample_output

  fam <- bed_fam_df(bo)
  ids <- bed_sample_IDs(bo)

  df <- fam[nrow(fam):1, ]
  annot <- fam[, 1:2]
  annot$SUBJID <- as.character(1:nrow(fam))
  annot$EXTRA <- 'useless'

  df2 <- annotate_plink_sample_output(df, annot, 'SUBJID')
  expect_equal(as.integer(df2$SUBJID), 89:1)

  expect_equivalent(df2[, -1], df[, -(1:2)])
}
test_that('annotate_plink_sample_output', .annotate_plink_sample_output())



.reorder_plink_snp_output <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  reorder_plink_snp_output <- plinker:::reorder_plink_snp_output

  bim <- bed_bim_df(bo)
  snp_ids <- bed_snp_IDs(bo)

  ### edge cases
  expect_error(reorder_plink_snp_output(iris, snp_ids), 'no SNP column')
  df <- bim
  df$SNP <- paste0('toto', df$SNP)
  expect_error(reorder_plink_snp_output(df, snp_ids), 'bad SNP ids')

  expect_identical(reorder_plink_snp_output(bim, NULL), bim)

  df2 <- reorder_plink_snp_output(bim, snp_ids)
  expect_identical(df2, bim)

  df <- bim[nrow(bim):1, ]
  df2 <- reorder_plink_snp_output(df, snp_ids)
  expect_identical(df2, bim)

  df3 <- reorder_plink_snp_output(bim, rev(snp_ids))
  expect_equivalent(df3, df)

  ### smaller sets
  df2 <- reorder_plink_snp_output(bim[1:10, ], rev(snp_ids))
  expect_equivalent(df2, bim[10:1, ])

  ### output may have several values by SNP (long format)
  df <- bim
  idx <- rep(1:nrow(bim), each = 4)
  df <- bim[idx, ]
  df$IDX <- 1:nrow(df)
  rownames(df) <- NULL

  df2 <- reorder_plink_snp_output(df, snp_ids)
  expect_identical(df2, df)

  df2 <- reorder_plink_snp_output(df, rev(snp_ids))
  expect_identical(unique(df2$SNP), rev(snp_ids))
}
test_that('reorder_plink_snp_output', .reorder_plink_snp_output())


