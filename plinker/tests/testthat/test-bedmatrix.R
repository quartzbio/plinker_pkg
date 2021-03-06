context('bedmatrix')

.bed_init_bedmatrix <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  mat <- plinker:::bed_init_bedmatrix(bo)

  expect_is(mat, 'BEDMatrix')

  expect_identical(rownames(mat), bed_sample_IDs(bo))
  expect_identical(colnames(mat), bed_snp_IDs(bo))
}
test_that('bed_init_bedmatrix', .bed_init_bedmatrix())



.bed_genotypes <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  ####################################
  mat <- bed_genotypes(bo)
  expect_is(mat, 'matrix')
  expect_true(is.integer(mat))
  expect_equal(dim(mat), c(bed_nb_samples(bo), bed_nb_snps(bo)))

  df <- bed_plink_freqx(bo, quiet = TRUE)
  df <- df[, c(2, 7:5, 10)]
  names(df) <- c("SNP", "A2A2", "A1A2", "A1A1", "NA")

  .mat_2_freq_df <- function(mat) {
    # trick: add one of each geno to get complete tables
    .ftable <- function(x) table(c(x, 0:2), useNA = "always")
    tbls <- apply(mat, 2, .ftable)
    df2 <- as.data.frame(t(tbls))
    df2 <- cbind(colnames(mat), df2, stringsAsFactors = FALSE)
    for (i in 2:4) df2[[i]] <- df2[[i]] - 1L
    names(df2) <- c("SNP", "A2A2", "A1A2", "A1A1", "NA")
    rownames(df2) <- NULL

    df2
  }

  df2 <- .mat_2_freq_df(mat)
  expect_identical(df2, df)

  ############################################
  mat2 <- bed_genotypes(bed_subset(bo, sample_idx = 11:20))
  expect_identical(mat2, mat[11:20, ])

  mat2 <- bed_genotypes(bed_subset(bo, snp_idx = 5:10))
  expect_identical(mat2, mat[, 5:10])

  mat2 <- bed_genotypes(bed_subset(bo, snp_idx = 5:10, sample_idx = 11:20))
  expect_identical(mat2, mat[11:20, 5:10])

  bo2 <- bed_subset(bed_subset(bo, snp_idx = 5:10, sample_idx = 11:20))
  expect_identical(bed_genotypes(bo2), mat[11:20, 5:10])

  expect_identical(bed_genotypes(bed_subset(bo, sample_idx = 1, snp_idx = 1)),
    mat[1, 1, drop = FALSE])

  ### ordering
  mat2 <- bed_genotypes(bed_subset(bo, snp_idx = 10:5, sample_idx = 20:11))
  expect_identical(mat2, mat[20:11, 10:5])

  snp_ids <- bed_snp_IDs(bo)
  sample_ids <- bed_sample_IDs(bo)
  mat2 <- bed_genotypes(bed_subset(bo, snp_IDs = snp_ids[c(5, 2)],
    sample_IDs = sample_ids[10:2]))
  expect_identical(mat2, mat[10:2, c(5, 2)])

  ### lexicographic_allele_order
  mat <- bed_genotypes(bo)
  mat2 <- bed_genotypes(bo, lexicographic_allele_order = TRUE)
  inv <- which(bed_allele2(bo) != bed_allele_higher(bo))
  expect_identical(mat2[, -inv], mat[, -inv])
  expect_identical(2L - mat2[, inv], mat[, inv])

  # with a subset
  bo2 <- bed_subset(bo, snp_idx = 5:10, sample_idx = 11:20)
  mat <- bed_genotypes(bo2)
  mat2 <- bed_genotypes(bo2, lexicographic_allele_order = TRUE)
  inv <- which(bed_allele2(bo2) != bed_allele_higher(bo2))
  expect_identical(mat2[, -inv], mat[, -inv])
  expect_identical(2L - mat2[, inv], mat[, inv])
}
test_that('bed_genotypes', .bed_genotypes())
