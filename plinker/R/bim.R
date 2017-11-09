#' read and parse a plink .bim file
#'
#' cf \url{http://www.cog-genomics.org/plink/1.9/formats#bim}
#'
#' @param path		the .bim file path
#' @return a data frame
#' @export
read_bim <- function(path) {
  df <- data.table::fread(path,
    sep = '\t',
    header = FALSE,
    verbose = FALSE,
    quote = '',
    strip.white = FALSE,
    data.table = FALSE,
    col.names = c('CHR', 'SNP', 'MORGANS', 'POS', 'A1', 'A2'),
    colClasses = c('character','character', 'numeric',
      'integer', 'character', 'character'),
    showProgress = FALSE
  )
}

#' save a plink .bim file
#'
#' @param bim_df		the .bim data frame
#' @param path			the .bim file path
#' @return a data frame
#' @export
#' @seealso read_bim
save_bim <- function(bim_df, path) {
  df <- readr::write_tsv(bim_df, path,
    col_names = FALSE)
}
