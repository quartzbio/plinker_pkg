context('genotypes')

.make_genotype_converters <- function() {
  make_genotype_converters <- plinker:::make_genotype_converters

  bo <- bed_open(plinker:::fetch_sample_bed())
  bim <- bed_bim_df(bo)

  convs <- make_genotype_converters('C', 'A')
  expect_equal(nrow(convs), 1)
  expect_equal(as.character(convs[1, 4:6, drop = TRUE]), c("A/A", "A/C", "C/C"))

  convs <- make_genotype_converters('C', 'A', sort = FALSE)
  expect_equal(as.character(convs[1, 4:6, drop = TRUE]), c("A/A", "C/A", "C/C"))

  convs <- make_genotype_converters('C', 'A', sep = '')
  expect_equal(as.character(convs[1, 4:6, drop = TRUE]), c("AA", "AC", "CC"))

  convs <- make_genotype_converters(bim$A1, bim$A2)
  expect_equal(dim(convs), c(11, 6))
}
test_that('make_genotype_converters', .make_genotype_converters())



.convert_genotypes_to_string <- function() {
  ### edge cases
  # convert_genotypes_to_string(c(0:3))


  genos <- matrix(c(0:2, NA, 2:0, NA), ncol = 1)

  expect_error(convert_genotypes_to_string(genos, LETTERS[1:2], LETTERS[1]),
    'bad param allele1')

  strs <- convert_genotypes_to_string(genos, 'C', 'A')
  expect_equal(dim(strs), dim(genos))
  expect_identical(strs[, 1], c("A/A", "A/C", "C/C", NA, "C/C", "A/C", "A/A", NA))

  bo <- bed_open(plinker:::fetch_sample_bed())
  genos <- bed_genotypes(bo)
  bim <- bed_bim_df(bo)

  ######################
  # N.B: sort=FALSE because plink does not sort het genotypes
  strs <- convert_genotypes_to_string(genos, bim$A1, bim$A2, sort = FALSE)
  expect_equal(dim(strs), dim(genos))
  expect_equal(dimnames(strs), dimnames(genos))

  setup_temp_dir()
  bed_plink_ped(bo, 'foo', quiet = TRUE)
  ped <- read_plink_ped('foo.ped')
  mat <- as.matrix(ped[, -(1:6)])
  # paste all pairs of columns

  .paste_pair <- function(i) {
    paste0(mat[, i], '/', mat[, i + 1])
  }
  strs2 <- sapply(seq.int(start = 1, length.out = ncol(mat)/2, by = 2), .paste_pair)
  strs2[strs2 == '0/0'] <- NA

  expect_equivalent(strs, strs2)
}
test_that('convert_genotypes_to_string', .convert_genotypes_to_string())

