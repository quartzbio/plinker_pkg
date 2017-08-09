
read_plink_output <- function(path, ...) {
  out <- utils::capture.output(
    df <- as.data.frame(
      readr::read_table(path, ..., progress = FALSE)
    )
  )

  df
}

read_plink_freq <- function(path) {
  read_plink_output(path, col_types = 'icccni')
}
