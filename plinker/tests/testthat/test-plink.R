context('plink')


.bed_plink_cmd <- function() {
  plink_cmd <- plinker:::plink_cmd

  bo <- bed_open(plinker:::fetch_sample_bed())

  setup_temp_dir()
  out <- bed_plink_cmd(bo, '--freq', quiet = TRUE)

  plink_frq <- 'plink.frq'

  expect_true(file.exists(plink_frq))
  df <- read_plink_freq(plink_frq)
  unlink(plink_frq)
  expect_equal(nrow(df), bed_nb_snps(bo))

  ### snp subset

  bed_plink_cmd(bo, snp_idx = 6:10, '--freq', quiet = TRUE)
  df2 <- read_plink_freq('plink.frq')
  bim_df <- bed_bim_df(bo, subset = FALSE)
  expect_identical(df2$SNP, bim_df[6:10, 'SNPID'])

  idx <- c(2, 10:16)
  bo2 <- bed_subset_snps_by_idx(bo, idx)

  bed_plink_cmd(bo2, '--freq', quiet = TRUE)
  df2 <- read_plink_freq('plink.frq')
  expect_identical(df2$SNP, bim_df[idx, 'SNPID'])

  ### sample subset
  sample_idx <- c(3, 11:20, 50:80)
  bed_plink_cmd(bo, sample_idx = sample_idx, '--freq', quiet = TRUE)
  df2 <- read_plink_freq('plink.frq')
  expect_equal(max(df2$NCHROBS), length(sample_idx) * 2)

  bo2 <- bed_subset_samples_by_idx(bo, sample_idx)
  bed_plink_cmd(bo2, '--freq', quiet = TRUE)
  df3 <- read_plink_freq('plink.frq')
  expect_identical(df3, df2)


  ### both subsets
  bo2 <- bed_subset_snps_by_idx(bo, idx)
  bo3 <- bed_subset_samples_by_idx(bo2, sample_idx)
  bed_plink_cmd(bo3, '--freq', quiet = TRUE)
  df4 <- read_plink_freq('plink.frq')

  expect_identical(df4$SNP, bim_df[idx, 'SNPID'])
  expect_equal(max(df4$NCHROBS), length(sample_idx) * 2)
}
test_that('bed_plink_cmd', .bed_plink_cmd())



.plink_cmd <- function() {
  plink_cmd <- plinker:::plink_cmd

  # edge cases
  suppressWarnings(
    expect_error(plink_cmd('toto', stderr = TRUE, stdout = TRUE),
      'ERROR running plink'))

  suppressWarnings(
    expect_error(plink_cmd('toto', command = '/not/found/foo',
         stderr = TRUE, stdout = TRUE),
      'error in running command'))

  text <- plink_cmd('--help', stdout = TRUE)
  expect_gt(length(text), 100)
}
test_that('plink_cmd', .plink_cmd())




.plink_version <- function() {
  ver <- plink_version()
  expect_match(ver, "PLINK v")
}
test_that('plink_version', .plink_version())

