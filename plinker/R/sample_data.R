get_extdata_dir <- function() {
  ### do not use devtools shimmed system.file
  dir <- base::system.file("extdata", package = 'plinker')
  if (identical(dir, '')) {
    dir <- base::system.file("inst", "extdata", package = 'plinker')
    stopifnot(nzchar(dir))
  }
  dir
}


fetch_sample_bed <- function() {
  file.path(get_extdata_dir(), 'plink')
}
