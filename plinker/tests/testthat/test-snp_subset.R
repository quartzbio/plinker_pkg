context('SNP subset')

.bed_subset_snps_by_idx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  nb <- bed_nb_snps(bo)
  expect_equal(nb, 17)

  ### edge cases
  expect_error(bed_subset_snps_by_idx(bo, NULL), 'empty snp_idx')
  expect_error(bed_subset_snps_by_idx(bo, 20), 'bad snp_idx range')

  bo2 <- bed_subset_snps_by_idx(bo, c(3, 5, 16, 3, 5))

  expect_equal(bed_nb_snps(bo2), 3)
  expect_equal(bed_nb_snps(bo2, subset = FALSE), 17)
  expect_identical(bed_snp_idx(bo2), c(3L, 5L, 16L))

  df <- bed_bim_df(bo2)
  expect_equal(rownames(df), as.character(c(3L, 5L, 16L)))
  expect_identical(df, bed_bim_df(bo2, subset = FALSE)[c(3L, 5L, 16L), ])

  ### recursive subsetting
  bo3 <- bed_subset_snps_by_idx(bo2, c(3, 1, 3))
  expect_equal(bed_nb_snps(bo3), 2)
  expect_equal(bed_nb_snps(bo3, subset = FALSE), 17)
  expect_identical(bed_snp_idx(bo3), c(3L, 16L))

  bo4 <-  bed_subset_snps_by_idx(bo3, 2)
  expect_equal(bed_nb_snps(bo4), 1)
  expect_identical(bed_snp_idx(bo4), 16L)
}
test_that('bed_subset_snps_by_idx', .bed_subset_snps_by_idx())



.bed_reset_subset_snps_by_idx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  expect_equal(bed_nb_snps(bo), 17)

  bo2 <- bed_subset_snps_by_idx(bo, c(3, 5, 16, 3, 5))
  expect_equal(bed_nb_snps(bo2), 3)

  bo3 <- bed_reset_subset_snps_by_idx(bo2)
  expect_equal(bed_nb_snps(bo3), 17)
}
test_that('bed_reset_subset_snps_by_idx', .bed_reset_subset_snps_by_idx())
