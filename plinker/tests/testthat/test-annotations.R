context('annotations')



.bed_set_snp_annot <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  expect_null(bed_get_snp_annot(bo))

  bim <- bed_bim_df(bo)
  ids <- bed_snp_IDs(bo)

  annot <- data.frame(
    SNP = rev(ids),
    ID = paste0('ID_', rev(ids)),
    POS = 1,
    stringsAsFactors = FALSE)

  expect_error(bed_set_snp_annot(bo, annot[, -1]), 'no SNP column')
  expect_error(bed_set_snp_annot(bo, annot[-1, ]), 'non exhaustive annotations')
  expect_error(bed_set_snp_annot(bo, rbind(annot, annot)),
     'duplicated values in SNP column')

  bo2 <- bed_set_snp_annot(bo, annot)

  annot2 <- bed_get_snp_annot(bo2)
  expect_equivalent(annot2[nrow(annot2):1, ], annot)


  ##### with secondary id -- -id=
  expect_error(bed_set_snp_annot(bo, annot, id = 'NOT'), 'bad column name')
  annot2 <- annot
  annot2$NUM <- nchar(annot2$ID)
  expect_error(bed_set_snp_annot(bo, annot2, id = 'NUM'), 'character')
  annot2$STR <- as.character(annot2$NUM)
  expect_warning(bed_set_snp_annot(bo, annot2, id = 'STR'), 'duplicated')

  bo2 <- bed_set_snp_annot(bo, annot, id = 'ID')
  bo3 <- bed_set_snp_annot(bo, annot)
  expect_identical(bed_get_snp_annot_id(bo2), 'ID')
  expect_null(bed_get_snp_annot_id(bo3))

  expect_identical(bed_get_snp_annot(bo3), bed_get_snp_annot(bo2))


  ### subsetting
  bo2 <- bed_set_snp_annot(bo, annot, 'ID')
  a2 <- bed_get_snp_annot(bo2)
  bo3 <- bed_subset(bo2, snp_idx = 10:5)
  a3 <- bed_get_snp_annot(bo3)

  expect_identical(a3, a2[10:5, ])
  expect_identical(bed_snp_IDs(bo3, custom = TRUE),
    bed_snp_IDs(bo2, custom = TRUE)[10:5])
}
test_that('bed_set_snp_annot', .bed_set_snp_annot())




.bed_set_get_sample_annot <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  expect_null(bed_get_sample_annot(bo))

  fam <- bed_fam_df(bo)
  ids <- bed_sample_IDs(bo)

  annot <- data.frame(
    MERGE_ID = rev(ids),
    SUBJID = paste0('ID_', rev(ids)),
    stringsAsFactors = FALSE)

  bo2 <- bed_set_sample_annot(bo, annot, rm_merge_id = FALSE)
  annot2 <- bed_get_sample_annot(bo2)
  ids2 <- bed_get_sample_annot(bo2, var = 'SUBJID')
  expect_identical(annot2$SUBJID, ids2)
  # MERGE_ID is still there
  expect_equivalent(annot2[nrow(annot2):1, ], annot)

  bo3 <- bed_set_sample_annot(bo, annot)
  annot3 <- bed_get_sample_annot(bo3)
  # no MERGE_ID anymore
  expect_identical(annot3, annot2[, -1, drop = FALSE])

  ids3 <- bed_get_sample_annot(bo2, var = 'SUBJID', fam = TRUE)
  # no OP
  expect_identical(ids3, ids2)
  annot3 <- bed_get_sample_annot(bo2, fam = TRUE)
  expect_identical(annot3[, 1:6], fam)


  ##### with secondary id -- -id=
  expect_error(bed_set_sample_annot(bo, annot, id = 'NOT'), 'bad column name')
  annot2 <- annot
  annot2$NUM <- nchar(annot2$SUBJID)
  expect_error(bed_set_sample_annot(bo, annot2, id = 'NUM'), 'character')
  annot2$STR <- as.character(annot2$NUM)
  expect_error(bed_set_sample_annot(bo, annot2, id = 'STR'), 'duplicated')

  bo2 <- bed_set_sample_annot(bo, annot, id = 'SUBJID')
  bo3 <- bed_set_sample_annot(bo, annot)
  expect_identical(bed_get_sample_annot_id(bo2), 'SUBJID')
  expect_null(bed_get_sample_annot_id(bo3))

  expect_identical(bed_get_sample_annot(bo3), bed_get_sample_annot(bo2))

  ### subsetting
  bo2 <- bed_set_sample_annot(bo, annot, id = 'SUBJID')
  a2 <- bed_get_sample_annot(bo2)
  bo3 <- bed_subset(bo2, sample_idx = 10:5)
  a3 <- bed_get_sample_annot(bo3)

  expect_identical(a3, a2[10:5, , drop = FALSE])
  expect_identical(bed_sample_IDs(bo3, custom = TRUE),
    bed_sample_IDs(bo2, custom = TRUE)[10:5])

}
test_that('bed_set_get_sample_annot', .bed_set_get_sample_annot())



