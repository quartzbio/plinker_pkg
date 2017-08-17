#' read a plink ped file
#'
#' @param path	the path of the ped file
#' @inheritParams params
#' @return the parsed ped as a data frame
#' @family read
#' @export
read_plink_ped <- function(path) {
  df <- read_plink_output(path)
  nsnps <- (ncol(df) - 6) / 2

  snp_cols <- rep(paste0('SNP', seq_len(nsnps)), each = 2)
  snp_cols <- paste0(snp_cols, c('A1', 'A2'))

  names(df) <- c('FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'SEX', 'PHENO', snp_cols)

  df
}

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
  df <- data.table::fread(path,
    verbose = FALSE,
    data.table = FALSE,
    showProgress = FALSE,
    ...
  )

  if (!is.null(df$CHR)) df$CHR <- as.character(df$CHR)

  df
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

