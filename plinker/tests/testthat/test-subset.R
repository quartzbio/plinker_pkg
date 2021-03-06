context('subsetting')



.bed_subset <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  snp_ids <- bed_bim_df(bo)$SNP
  sample_ids <- bed_sample_IDs(bo)

  ### edge cases
  expect_error(bed_subset(bo, snp_IDs = snp_ids, snp_idx = 1),
    "incompatible parameters")
  expect_error(bed_subset(bo, sample_IDs = sample_ids, sample_idx = 1),
    "incompatible parameters")
  expect_error(bed_subset(bo, sample_idx = 1, snp_idx = 1, invert = TRUE),
    "can not use both SNP and sample subsetting with invert=TRUE")


  expect_identical(bed_subset(bo, snp_IDs = snp_ids[10:15]),
    bed_subset_snps_by_IDs(bo, snp_ids[10:15]))
  expect_identical(bed_subset(bo, snp_idx = 3),
    bed_subset_snps_by_idx(bo, 3))
  expect_identical(bed_subset(bo, sample_IDs = sample_ids[10:15]),
    bed_subset_samples_by_IDs(bo, sample_ids[10:15]))
  expect_identical(bed_subset(bo, sample_idx = 3),
    bed_subset_samples_by_idx(bo, 3))

  # invert
  expect_identical(bed_subset(bo, snp_IDs = snp_ids[1:10], invert = TRUE),
    bed_subset_snps_by_IDs(bo, snp_ids[11:17]))
  expect_identical(bed_subset(bo, snp_idx = 17:11, invert = TRUE),
    bed_subset_snps_by_idx(bo, 1:10))
  expect_identical(bed_subset(bo, sample_IDs = sample_ids[1:10], invert = TRUE),
    bed_subset_samples_by_IDs(bo, sample_ids[11:89]))
  expect_identical(bed_subset(bo, sample_idx = 81:89, invert = TRUE),
    bed_subset_samples_by_idx(bo, 1:80))

  bo2 <- bed_subset(bo, snp_idx = 3:9, sample_IDs = sample_ids[10:15])
  expect_identical(bed_snp_idx(bo2), 3:9)
  expect_identical(bed_sample_idx(bo2), 10:15)

  bo3 <- bed_subset(bo2, snp_IDs = snp_ids[5:6], sample_idx = 4:2)
  expect_identical(bed_snp_idx(bo3), 5:6)
  expect_identical(bed_sample_idx(bo3), (10:15)[4:2])

  # invert
  bo3 <- bed_subset(bo2, snp_IDs = snp_ids[5:6], invert = TRUE)
  expect_identical(bed_snp_IDs(bo3), snp_ids[c(3:4, 7:9)])
  bo3 <- bed_subset(bo2, sample_idx = 4:2, invert = TRUE)
  expect_identical(bed_sample_idx(bo3), c(10L, 14L, 15L))
}
test_that('bed_subset', .bed_subset())



.bed_subset_samples_by_idx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  nb <- bed_nb_samples(bo)
  expect_equal(nb, 89)

  ### edge cases
  expect_error(bed_subset_samples_by_idx(bo, NULL), 'empty sample_idx')
  expect_error(bed_subset_samples_by_idx(bo, 100), 'bad sample_idx range')
  # duplicates
  expect_error(bed_subset_samples_by_idx(bo, c(1, 2, 1)), 'duplicated')
  expect_error(bed_subset_samples_by_idx(bo, c(1, NA)), 'missing')

  bo2 <- bed_subset_samples_by_idx(bo, c(21:35))

  expect_equal(bed_nb_samples(bo2), 15)
  expect_equal(bed_nb_samples(bo2, subset = FALSE), 89)
  expect_identical(bed_sample_idx(bo2), 21:35)

  df <- bed_fam_df(bo2)
  expect_equal(rownames(df), as.character(21:35))
  expect_identical(df, bed_fam_df(bo2, subset = FALSE)[21:35, ])

  # invert N.B: ordering is meaningless
  bo2 <- bed_subset_samples_by_idx(bo, 89:11, invert = TRUE)
  expect_identical(bed_sample_idx(bo2), 1:10)

  ### recursive subsetting
  bo2 <- bed_subset_samples_by_idx(bo, c(21:35))
  bo3 <- bed_subset_samples_by_idx(bo2, c(3, 1, 5))
  expect_equal(bed_nb_samples(bo3), 3)
  expect_equal(bed_nb_samples(bo3, subset = FALSE), 89)
  expect_identical(bed_sample_idx(bo3), c(23L, 21L, 25L))

  bo4 <-  bed_subset_samples_by_idx(bo3, 2)
  expect_equal(bed_nb_samples(bo4), 1)
  expect_identical(bed_sample_idx(bo4), 21L)

  # invert
  bo4 <- bed_subset_samples_by_idx(bo3, 2:3, invert = TRUE)
  expect_identical(bed_sample_idx(bo4), 23L)

  ### ordering
  bo2 <- bed_subset_samples_by_idx(bo, 10:5)
  expect_identical(bed_sample_idx(bo2), 10:5)
}
test_that('bed_subset_samples_by_idx', .bed_subset_samples_by_idx())



.bed_reset_subset_samples_by_idx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  expect_equal(bed_nb_samples(bo), 89)

  bo2 <- bed_subset_samples_by_idx(bo, c(3, 5, 16))
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

  expect_error(bed_sample_IDs_to_idx(bo, NULL), 'empty sample IDs')
  expect_error(bed_sample_IDs_to_idx(bo, sample_ids[c(1,1)]), 'duplicated')

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
  expect_error(bed_subset_samples_by_IDs(bo, NULL), 'empty sample IDs')
  expect_error(bed_subset_samples_by_IDs(bo, "toto"), 'bad sample ID')

  ids <- bed_sample_IDs(bo)

  bo2 <- bed_subset_samples_by_IDs(bo, ids[40:30])

  expect_equal(bed_nb_samples(bo2), 11)
  expect_equal(bed_nb_samples(bo2, subset = FALSE), 89)
  expect_identical(bed_sample_idx(bo2), 40:30)

  ### recursive subsetting
  ids3 <- ids[35:38]
  bo3 <- bed_subset_samples_by_IDs(bo2, ids3)
  expect_equal(bed_nb_samples(bo3), 4)
  expect_equal(bed_nb_samples(bo3, subset = FALSE), 89)
  expect_identical(bed_sample_idx(bo3), 35:38)
  expect_identical(bed_sample_IDs(bo3), ids3)

  bo4 <- bed_subset_samples_by_IDs(bo3, ids3[2])
  expect_equal(bed_nb_samples(bo4), 1)
  expect_identical(bed_sample_idx(bo4), 36L)
  expect_identical(bed_sample_IDs(bo4), ids3[2])

  bo4 <- bed_subset_samples_by_IDs(bo3, ids3[1:2], invert = TRUE)
  expect_identical(bed_sample_IDs(bo4), ids3[3:4])
}
test_that('bed_subset_samples_by_IDs', .bed_subset_samples_by_IDs())


.bed_subset_snps_by_IDs <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  nb <- bed_nb_snps(bo)
  expect_equal(nb, 17)

  ### edge cases
  expect_error(bed_subset_snps_by_IDs(bo, NULL), 'empty snp_idx')
  expect_error(bed_subset_snps_by_IDs(bo, "toto"), 'bad snp_IDs')

  bim_df <- bed_bim_df(bo)
  idx <- c(5:12, 1L, 16L)
  ids <- bim_df$SNP[idx]

  bo2 <- bed_subset_snps_by_IDs(bo, ids)

  expect_equal(bed_nb_snps(bo2), 10)
  expect_equal(bed_nb_snps(bo2, subset = FALSE), 17)

  expect_identical(bed_snp_idx(bo2), idx)
  expect_identical(bed_snp_IDs(bo2), ids)

  df <- bed_bim_df(bo2)
  expect_equal(rownames(df), as.character(idx))
  expect_identical(df, bed_bim_df(bo2, subset = FALSE)[idx, ])

  # invert
  bo2 <- bed_subset_snps_by_IDs(bo, ids, invert = TRUE)
  expect_identical(bed_snp_IDs(bo2), setdiff(bed_snp_IDs(bo), ids))

  ### recursive subsetting
  ids2 <- ids[c(2, 4:7, 9)]
  bo2 <- bed_subset_snps_by_IDs(bo, ids)
  bo3 <- bed_subset_snps_by_IDs(bo2, ids2)
  expect_equal(bed_nb_snps(bo3), length(ids2))
  expect_equal(bed_nb_snps(bo3, subset = FALSE), 17)
  expect_identical(sort(bim_df$SNP[bed_snp_idx(bo3)]), sort(ids2))
  expect_identical(sort(bed_snp_IDs(bo3)), sort(unique(ids2)))

  bo4 <- bed_subset_snps_by_IDs(bo3, "rs10105623")
  expect_equal(bed_nb_snps(bo4), 1)
  expect_identical(bim_df$SNP[bed_snp_idx(bo4)], "rs10105623")
  expect_identical(bed_snp_IDs(bo4), "rs10105623")

  expect_identical(bed_snp_IDs(bo4, subset = FALSE), bim_df$SNP)

  # invert
  bo4 <- bed_subset_snps_by_IDs(bo3, c("rs7835221", "rs2460911", "rs12156420"),
    invert = TRUE)
  expect_identical(bed_snp_IDs(bo4), c("rs10105623", "rs17786052", "rs17121574"))
}
test_that('bed_subset_snps_by_IDs', .bed_subset_snps_by_IDs())



.bed_subset_snps_by_idx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  nb <- bed_nb_snps(bo)
  expect_equal(nb, 17)

  ### edge cases
  expect_error(bed_subset_snps_by_idx(bo, NULL), 'empty snp_idx')
  expect_error(bed_subset_snps_by_idx(bo, 20), 'bad snp_idx range')
  expect_error(bed_subset_snps_by_idx(bo, c(1, 2, 1)), 'duplicated')
  expect_error(bed_subset_snps_by_idx(bo, c(1, NA)), 'missing')

  # invert
  bo2 <- bed_subset_snps_by_idx(bo, 11:17, invert = TRUE)
  expect_identical(bed_snp_idx(bo2), 1:10)

  bo2 <- bed_subset_snps_by_idx(bo, c(3, 5, 16))

  expect_equal(bed_nb_snps(bo2), 3)
  expect_equal(bed_nb_snps(bo2, subset = FALSE), 17)
  expect_identical(bed_snp_idx(bo2), c(3L, 5L, 16L))

  df <- bed_bim_df(bo2)
  expect_equal(rownames(df), as.character(c(3L, 5L, 16L)))
  expect_identical(df, bed_bim_df(bo2, subset = FALSE)[c(3L, 5L, 16L), ])


  ### recursive subsetting
  bo2 <- bed_subset_snps_by_idx(bo, c(3, 5, 16))
  bo3 <- bed_subset_snps_by_idx(bo2, c(3, 1))
  expect_equal(bed_nb_snps(bo3), 2)
  expect_equal(bed_nb_snps(bo3, subset = FALSE), 17)
  expect_identical(bed_snp_idx(bo3), c(16L, 3L))

  bo4 <- bed_subset_snps_by_idx(bo3, 2)
  expect_equal(bed_nb_snps(bo4), 1)
  expect_identical(bed_snp_idx(bo4), 3L)

  # invert
  bo4 <- bed_subset_snps_by_idx(bo3, 1, invert = TRUE)
  expect_identical(bed_snp_idx(bo4), 3L)
}
test_that('bed_subset_snps_by_idx', .bed_subset_snps_by_idx())



.bed_reset_subset_snps_by_idx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  expect_equal(bed_nb_snps(bo), 17)

  bo2 <- bed_subset_snps_by_idx(bo, c(3, 5, 16))
  expect_equal(bed_nb_snps(bo2), 3)

  bo3 <- bed_reset_subset_snps_by_idx(bo2)
  expect_equal(bed_nb_snps(bo3), 17)
}
test_that('bed_reset_subset_snps_by_idx', .bed_reset_subset_snps_by_idx())
