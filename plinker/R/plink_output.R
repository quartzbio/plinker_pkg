
read_plink_output <- function(path, ...) {
#  df <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
#  if (!is.null(df$CHR)) df$CHR <- as.character(df$CHR)
#
#  df
#  out <- utils::capture.output(
#    df <- as.data.frame(
#      readr::read_table(path, ..., progress = FALSE)
#    )
#  )
#
#  df
  data.table::fread(path,
    verbose = FALSE,
    data.table = FALSE,
    showProgress = FALSE,
    ...
  )
}

read_plink_freq <- function(path) {
  read_plink_output(path)
  #read_plink_output(path, col_types = 'ccccni')
}

read_plink_freq_counts <- function(path) {
  read_plink_output(path,
    colClasses = c('character','character', 'character', 'character',
      'integer', 'integer', 'integer'))
#  read_plink_output(path, col_types = 'cccciii')
}

