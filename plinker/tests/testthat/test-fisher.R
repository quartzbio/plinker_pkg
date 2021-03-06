context('fisher test')


.bed_plink_fisher <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)
  fam <- bed_fam_df(bo)

  ### edge cases
  expect_error(bed_plink_fisher(bo, phenotype = 1L, quiet = TRUE),
    'bad phenotype length')
  expect_error(bed_plink_fisher(bo, phenotype = as.numeric(fam$PHENO),
      quiet = TRUE), 'must be an integer vector')
  pheno <- fam$PHENO
  pheno[5] <- NA
  # NA are automatically recoded
  expect_error(bed_plink_fisher(bo, phenotype = pheno, quiet = TRUE), NA)


  expect_error(bed_plink_fisher(bo, phenotype = rep(3L, length(pheno)),
       quiet = TRUE), 'no non-missing phenotypes')

  out <- bed_plink_fisher(bo, quiet = TRUE)

  expect_true(nrow(out) == nsnp * 5) # 5 tests are performed for each snp
  expect_identical(colnames(out), c("CHR", "SNP", "A1", "A2", "TEST",
          "AFF", "UNAFF", "P"))

   ### only check for the GENO TEST
  .check_on_geno <- function(out, ref) {
    outgeno <- out[out$TEST == 'GENO', ]
    rownames(outgeno) <- NULL
    expect_equal(outgeno, ref, tolerance = 1e-3)
  }

  .check_on_geno(out, bed_R_fisher(bo))

  out2 <- bed_plink_fisher(bo, phenotype = fam$PHENO, quiet = TRUE)
  expect_identical(out2, out)

  out2 <- bed_plink_fisher(bo, phenotype = rev(fam$PHENO), quiet = TRUE)
  expect_false(identical(out2, out))

  ### missing_phenotype
  pheno <- fam$PHENO
  nmiss <- function(aff, unaff) {
    sum(as.integer(unlist(strsplit(c(aff, unaff), '/'))))
  }

  out <- bed_plink_fisher(bo, pheno, quiet = TRUE)
  expect_equal(nmiss(out$AFF[1], out$UNAFF[1]), 89)

  pheno[11:20] <- NA
  out2 <- bed_plink_fisher(bo, pheno, quiet = TRUE)
  expect_equal(nmiss(out2$AFF[1], out2$UNAFF[1]), 79)

  pheno[11:20] <- -9L # plink encoding value
  # pheno has 3 values since the -9 is not identified as the missing phenotype
  expect_error(bed_plink_fisher(bo, pheno, quiet = TRUE),
    'phenotype must be binary')

  ### ordering
  bo2 <- bed_subset(bo, snp_idx = 10:2)
  out <- bed_plink_fisher(bo2, quiet = TRUE)
  expect_identical(unique(out$SNP), bed_snp_IDs(bo2))

  ### allele order: inversed OR, opposite stat, same pvalue
  ### PB: bug in PLINK, lots of NA
  res1 <- bed_plink_fisher(bo, quiet = TRUE)
  res2 <- bed_plink_fisher(bo, quiet = TRUE, lexicographic_allele_order = TRUE)
  res3 <- bed_R_fisher(bo, lexicographic_allele_order = TRUE)

  res2geno <- res2[res2$TEST == 'GENO', ]
  rownames(res2geno) <- NULL
  expect_equal(res2geno, res3, tolerance = 1e-4)

  # non swapped
  expect_identical(res2[1:5, ], res1[1:5, ])
  # swapped
  df1 <- res1[6:10, ]
  df2 <- res2[6:10, ]
  expect_false(any(df1$AFF == df2$AFF))
  expect_false(any(df1$UNAFF == df2$UNAFF))
  expect_equal(df2$P[1:3], df1$P[1:3], tolerance = 1e-4)
  # recessive and dominant are reversed
  expect_equal(df2$P[4:5], df1$P[5:4], tolerance = 1e-4)

  ### snp annot
  bim <- bed_bim_df(bo)
  annot <- data.frame(
    SNP = bim$SNP,
    ID = paste0('ID_', bim$SNP),
    TOTO = 1,
    stringsAsFactors = FALSE)

  bo2 <- bed_set_snp_annot(bo, annot, 'ID')
  res <- bed_plink_fisher(bo, quiet = TRUE)
  res2 <- bed_plink_fisher(bo2, quiet = TRUE)

  expect_identical(res2[, -3], res)
  expect_identical(paste0('ID_', res2$SNP), res2$ID)
}
test_that('bed_plink_fisher', .bed_plink_fisher())

