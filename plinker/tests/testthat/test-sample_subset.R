context('Sample subset')

.bed_subset_samples_by_idx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  nb <- bed_nb_samples(bo)
  expect_equal(nb, 89)

  ### edge cases
  expect_error(bed_subset_samples_by_idx(bo, NULL), 'empty sample_idx')
  expect_error(bed_subset_samples_by_idx(bo, 100), 'bad sample_idx range')

  bo2 <- bed_subset_samples_by_idx(bo, c(21:30, 25:35))

  expect_equal(bed_nb_samples(bo2), 15)
  expect_equal(bed_nb_samples(bo2, subset = FALSE), 89)
  expect_identical(bed_sample_idx(bo2), 21:35)

  df <- bed_fam_df(bo2)
  expect_equal(rownames(df), as.character(21:35))
  expect_identical(df, bed_fam_df(bo2, subset = FALSE)[21:35, ])

  ### recursive subsetting
  bo3 <- bed_subset_samples_by_idx(bo2, c(3, 1, 3))
  expect_equal(bed_nb_samples(bo3), 2)
  expect_equal(bed_nb_samples(bo3, subset = FALSE), 89)
  expect_identical(bed_sample_idx(bo3), c(21L, 23L))

  bo4 <-  bed_subset_samples_by_idx(bo3, 2)
  expect_equal(bed_nb_samples(bo4), 1)
  expect_identical(bed_sample_idx(bo4), 23L)
}
test_that('bed_subset_samples_by_idx', .bed_subset_samples_by_idx())



.bed_reset_subset_samples_by_idx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  expect_equal(bed_nb_samples(bo), 89)

  bo2 <- bed_subset_samples_by_idx(bo, c(3, 5, 16, 3, 5))
  expect_equal(bed_nb_samples(bo2), 3)

  bo3 <- bed_reset_subset_samples_by_idx(bo2)
  expect_equal(bed_nb_samples(bo3), 89)
}
test_that('bed_reset_subset_samples_by_idx', .bed_reset_subset_samples_by_idx())
