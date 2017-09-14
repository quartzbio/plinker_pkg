context('plink')


.save_plink_with_lexicographic_alleles_order <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  a1 <- bed_allele1(bo)
  a2 <- bed_allele2(bo)

  # not ordered
  expect_equal(sum(a1 > a2), 9)

  setup_temp_dir()
  save_plink_with_lexicographic_alleles_order(plinker:::fetch_sample_bed(),
    'toto', quiet = F)

  bo2 <- bed_open('toto')
  b1 <- bed_allele1(bo2)
  b2 <- bed_allele2(bo2)
  # check alleles
  expect_equal(sum(b1 > b2), 0)

  expect_identical(bed_genotypes_as_strings(bo2), bed_genotypes_as_strings(bo))
}
test_that('save_plink_with_lexicographic_alleles_order',
  .save_plink_with_lexicographic_alleles_order())



.save_plink_alleles <- function() {
  save_plink_alleles <- plinker:::save_plink_alleles
  bo <- bed_open(plinker:::fetch_sample_bed())

  bim <- bed_bim_df(bo)

  df <- bim[, c('SNP', 'A1', 'A2')]
  setup_temp_dir()

  save_plink_alleles(df, 'toto.tall')
  expect_true(file.exists('toto.tall'))

  df2 <- plinker:::read_plink_output('toto.tall', header = FALSE)
  names(df2) <- names(df)
  expect_identical(df2, df)
}
test_that('save_plink_alleles', .save_plink_alleles())



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



.plink_bfile_to_tfile_tobfile <- function() {
  bfile_prefix <- plinker:::fetch_sample_bed()
  setup_temp_dir()
  bo1 <- bed_open(bfile_prefix)

  tfile_prefix <- 'plinkt'
  plink_bfile_to_tfile(bfile_prefix, tfile_prefix, quiet = TRUE)
  expect_true(all(file.exists(paste0('plinkt.', c('tped', 'tfam')))))


  plink_tfile_to_bfile(tfile_prefix, 'toto', quiet = TRUE)
  expect_true(all(file.exists(paste0('toto.', c('bed', 'fam', 'bim')))))

  # compare them
  bo1 <- bed_open(bfile_prefix)
  bo2 <- bed_open('toto')

  g1 <- bed_genotypes_as_strings(bo1)
  g2 <- bed_genotypes_as_strings(bo2)
  expect_identical(g1, g2)

  df1 <- bed_bim_df(bo1)
  df2 <- bed_bim_df(bo2)

  swaps <- which(df1$A1 != df2$A1)
  expect_equal(swaps, 15)
  expect_identical(df1$A1[swaps], df2$A2[swaps])
  expect_identical(df1$A2[swaps], df2$A1[swaps])

  df1_fixed <- df1
  df1_fixed$A1[swaps] <- df2$A1[swaps]
  df1_fixed$A2[swaps] <- df2$A2[swaps]
  expect_identical(df1_fixed, df2)

  g1 <- bed_genotypes(bo1)
  g2 <- bed_genotypes(bo2)

  g1[, swaps] <- 2L - g1[, swaps]
  expect_identical(g1, g2)
}
test_that('plink_bfile_to_tfile_tobfile', .plink_bfile_to_tfile_tobfile())