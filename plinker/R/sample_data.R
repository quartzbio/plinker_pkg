get_extdata_dir <- function(debug = FALSE) {
  ### do not use devtools shimmed system.file
  dir <- base::system.file("extdata", package = 'plinker')
  if (debug || identical(dir, '')) {
    dir <- base::system.file("inst", "extdata", package = 'plinker')
    stopifnot(nzchar(dir))
  }
  dir
}


fetch_sample_bed <- function() {
  file.path(get_extdata_dir(), 'plink')
}
