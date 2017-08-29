context('R stats')

.bed_R_lm <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  bo1 <- bed_subset(bo, snp_idx = 10)
  res <- bed_R_lm(bo1)


  browser()
}
test_that('bed_R_lm', .bed_R_lm())

