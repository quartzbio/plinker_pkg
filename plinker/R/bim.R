#' read and parse a plink .bim file
#'
#' cf \url{http://www.cog-genomics.org/plink/1.9/formats#bim}
#'
#' @param path		the .bim file path
#' @return a data frame
#' @export
read_bim <- function(path) {
  as.data.frame(readr::read_delim(path, '\t',
      col_names = c('CHR', 'SNPID', 'MORGANS', 'POS', 'A1', 'A2'),
      col_types = 'ccnicc',
      progress = FALSE
   ))
}

