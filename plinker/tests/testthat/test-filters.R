context('filters')

.bed_filter_snps_by_missing_rate <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  lmiss <- bed_plink_missing(bo, quiet = TRUE)$snp
  expect_false(all(lmiss$F_MISS < 2/89))

  bo2 <- bed_filter_snps_by_missing_rate(bo, 2/89, quiet = TRUE)
  expect_equal(bed_nb_snps(bo2), bed_nb_snps(bo) - 1)

  lmiss <- bed_plink_missing(bo2, quiet = TRUE)$snp
  expect_true(all(lmiss$F_MISS < 2/89))
}
test_that('bed_filter_snps_by_missing_rate', .bed_filter_snps_by_missing_rate())

