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



.bed_sample_IDs_to_idx <- function() {
  bed_sample_IDs_to_idx <- plinker:::bed_sample_IDs_to_idx

  bo <- bed_open(plinker:::fetch_sample_bed())
  fam_df <- bed_fam_df(bo)

  sample_ids <- bed_sample_IDs(bo)
  sample_iids <- bed_sample_IDs(bo, ignore_fid = TRUE)

  expect_error(bed_sample_IDs_to_idx(bo, fam_df$IID), 'bad sample IDs')

  expect_identical(bed_sample_IDs_to_idx(bo, rev(sample_ids)), nrow(fam_df):1)
  expect_identical(bed_sample_IDs_to_idx(bo, sample_iids, ignore_fid = TRUE),
    1:nrow(fam_df))

  expect_identical(bed_sample_IDs_to_idx(bo, sample_ids[51]), 51L)


  ### subsetting
  idx <- c(21:40, 7, 69)
  bo2 <- bed_subset_samples_by_idx(bo, idx)

  # out of range
  expect_error(bed_sample_IDs_to_idx(bo2, sample_ids), 'bad sample IDs')

  expect_identical(bed_sample_IDs_to_idx(bo2, bed_sample_IDs(bo2)[13:5]), 13:5)

  expect_error(bed_sample_IDs_to_idx(bo2,
      bed_sample_IDs(bo2, ignore_fid = TRUE)[13:5]), 'bad sample IDs')
  expect_identical(bed_sample_IDs_to_idx(bo2,
      bed_sample_IDs(bo2, ignore_fid = TRUE)[13:5], ignore_fid = TRUE), 13:5)
}
test_that('bed_sample_IDs_to_idx', .bed_sample_IDs_to_idx())



.bed_subset_samples_by_IDs <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  nb <- bed_nb_samples(bo)
  expect_equal(nb, 89)

  ### edge cases
  expect_error(bed_subset_samples_by_IDs(bo, NULL), 'empty sample_idx')
  expect_error(bed_subset_samples_by_IDs(bo, "toto"), 'bad sample ID')

  ids <- bed_sample_IDs(bo)

  bo2 <- bed_subset_samples_by_IDs(bo, ids[40:30])

  expect_equal(bed_nb_samples(bo2), 11)
  expect_equal(bed_nb_samples(bo2, subset = FALSE), 89)
  expect_identical(bed_sample_idx(bo2), 30:40)

  ### recursive subsetting
  ids3 <- ids[35:38]
  bo3 <- bed_subset_samples_by_IDs(bo2, ids3)
  expect_equal(bed_nb_samples(bo3), 4)
  expect_equal(bed_nb_samples(bo3, subset = FALSE), 89)
  expect_identical(bed_sample_idx(bo3), 35:38)
  expect_identical(bed_sample_IDs(bo3), ids3)

  bo4 <-  bed_subset_samples_by_IDs(bo3, ids3[2])
  expect_equal(bed_nb_samples(bo4), 1)
  expect_identical(bed_sample_idx(bo4), 36L)
  expect_identical(bed_sample_IDs(bo4), ids3[2])
}
test_that('bed_subset_samples_by_IDs', .bed_subset_samples_by_IDs())
