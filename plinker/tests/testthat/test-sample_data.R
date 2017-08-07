context('sample_data')

.fetch_sample_bed <- function() {
  ds <- plinker:::fetch_sample_bed()
  expect_true(all(file.exists(plinker:::unprefix_bed(ds))))
}
test_that('fetch_sample_bed', .fetch_sample_bed())

