context('plink')


.bed_plink_missing <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  fam <- bed_fam_df(bo)
  bim <- bed_bim_df(bo)

  res <- bed_plink_missing(bo, quiet = TRUE)

  imiss <- res$sample
  expect_identical(imiss[, 1:2], fam[, 1:2])
  expect_equal(unique(imiss$N_GENO), 17)
  expect_equal(max(imiss$N_MISS), 1)
  expect_equal(max(imiss$F_MISS), 1/17, tolerance = 1e-5)

  lmiss <- res$snp

  expect_identical(lmiss[, 1:2], bim[, 1:2])
  gfreqs <- bed_plink_freqx(bo, quiet = TRUE)
  expect_equal(gfreqs[[10]], lmiss$N_MISS)
  expect_equal(gfreqs[[10]] / 89, lmiss$F_MISS, tolerance = 1e-5)
}
test_that('bed_plink_missing', .bed_plink_missing())




.bed_plink_freqx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  bim_df <- bed_bim_df(bo)

  ### no subset
  df <- bed_plink_freqx(bo, quiet = TRUE)

  nb <- unique(rowSums(df[, 5:10]))
  expect_equal(nb, bed_nb_samples(bo))

  expect_identical(df$SNP, bim_df$SNP)
  expect_identical(df$A1, bim_df$A1)
  expect_identical(df$A2, bim_df$A2)

  # sample filtering
  df2 <- bed_plink_freqx(bo, allow_no_sex = FALSE, quiet = TRUE)
  # only ignore phenotype for no sex, does not impact the freqs
  expect_identical(df2, df)

  df2 <- bed_plink_freqx(bo, nonfounders = FALSE, quiet = TRUE)
  nb <- unique(rowSums(df2[, 5:10]))
  nb_nonfounders <- bed_nb_samples(bo) - nb
  expect_equal(nb_nonfounders, 3)

  ## !! 3 non founders, although only 2 actual ones since on hasa parent not
  ## in the fam file
}
test_that('bed_plink_freqx', .bed_plink_freqx())



.bed_plink_freq_count <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  bim_df <- bed_bim_df(bo)

  ### no subset
  df <- bed_plink_freq_count(bo, quiet = TRUE)

  nb <- unique(with(df, C1 + C2 + 2*G0)) / 2
  expect_equal(nb, bed_nb_samples(bo))

  expect_identical(df$SNP, bim_df$SNP)
  expect_identical(df$A1, bim_df$A1)
  expect_identical(df$A2, bim_df$A2)

  # sample filtering
  df2 <- bed_plink_freq_count(bo, allow_no_sex = FALSE, quiet = TRUE)
  # only ignore phenotype for no sex, does not impact the freqs
  expect_identical(df2, df)

  df2 <- bed_plink_freq_count(bo, nonfounders = FALSE, quiet = TRUE)
  nb <- unique(with(df2, C1 + C2 + 2*G0)) / 2
  nb_nonfounders <- bed_nb_samples(bo) - nb
  expect_equal(nb_nonfounders, 3)

  ## !! 3 non founders, although only 2 actual ones since on hasa parent not
  ## in the fam file

  #### subset by snps
  bo2 <- bed_subset_snps_by_idx(bo, 6:10)
  df2 <- bed_plink_freq_count(bo2, quiet = TRUE)
  expect_equivalent(df2, df[6:10, ])

  ### subset by samples
  bo2 <- bed_subset_samples_by_idx(bo, 21:60)
  df2 <- bed_plink_freq_count(bo2, quiet = TRUE)
  nb <- unique(with(df2, C1 + C2 + 2*G0)) / 2
  expect_equal(nb, 40)
}
test_that('bed_plink_freq_count', .bed_plink_freq_count())



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
  expect_identical(df2$SNP, bim_df[6:10, 'SNP'])

  idx <- c(2, 10:16)
  bo2 <- bed_subset_snps_by_idx(bo, idx)

  bed_plink_cmd(bo2, '--freq', quiet = TRUE)
  df2 <- read_plink_freq('plink.frq')
  expect_identical(df2$SNP, bim_df[idx, 'SNP'])

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

  expect_identical(df4$SNP, bim_df[idx, 'SNP'])
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

