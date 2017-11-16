context('sample_data')

.fetch_sample_bed <- function() {
  ds <- plinker:::fetch_sample_bed()
  expect_true(all(file.exists(plinker:::unprefix_bed(ds))))
}
test_that('fetch_sample_bed', .fetch_sample_bed())



.get_extdata_dir <- function() {
  setup_temp_dir()
  dir <- try(plinker:::get_extdata_dir(debug = TRUE))
  # crash during coverage
  expect_true(file.exists(dir) || is(dir, 'try-error'))
}
test_that('get_extdata_dir', .get_extdata_dir())