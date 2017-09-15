context('plink output')


.reorder_plink_output <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  reorder_plink_output <- plinker:::reorder_plink_output

  bim <- bed_bim_df(bo)
  snp_ids <- bed_snp_IDs(bo)

  ### edge cases
  expect_error(reorder_plink_output(iris, snp_ids), 'no SNP column')
  df <- bim
  df$SNP <- paste0('toto', df$SNP)
  expect_error(reorder_plink_output(df, snp_ids), 'bad SNP ids')

  expect_identical(reorder_plink_output(bim, NULL), bim)

  df2 <- reorder_plink_output(bim, snp_ids)
  expect_identical(df2, bim)

  df <- bim[nrow(bim):1, ]
  df2 <- reorder_plink_output(df, snp_ids)
  expect_identical(df2, bim)

  df3 <- reorder_plink_output(bim, rev(snp_ids))
  expect_equivalent(df3, df)

  ### smaller sets
  df2 <- reorder_plink_output(bim[1:10, ], rev(snp_ids))
  expect_equivalent(df2, bim[10:1, ])

  ### output may have several values by SNP (long format)
  df <- bim
  idx <- rep(1:nrow(bim), each = 4)
  df <- bim[idx, ]
  df$IDX <- 1:nrow(df)
  rownames(df) <- NULL

  df2 <- reorder_plink_output(df, snp_ids)
  browser()
  expect_identical(df2, df)

  df2 <- reorder_plink_output(df, rev(snp_ids))
  expect_identical(unique(df2$SNP), rev(snp_ids))




}
test_that('reorder_plink_output', .reorder_plink_output())


