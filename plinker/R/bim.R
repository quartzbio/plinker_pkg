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
    col.names = c('CHR', 'SNPID', 'MORGANS', 'POS', 'A1', 'A2'),
    colClasses = c('character','character', 'numeric',
      'integer', 'character', 'character'),
    showProgress = FALSE
  )
#  as.data.frame(readr::read_delim(path, '\t',
#      col_names = c('CHR', 'SNPID', 'MORGANS', 'POS', 'A1', 'A2'),
#      col_types = 'ccnicc',
#      progress = FALSE
#   ))
}


#read_bim_snpids_by_block <- function(path, blocks) {
#
#  nb <- nrow(blocks)
#  starts <- blocks[, 1]
#  widths <- blocks[, 2]
#  ends <- starts + widths - 1L
#  skips <- starts - c(0, ends[-nb]) - 1L
#
#  con <- file(path, 'rb')
#  open(con)
#  read_delim <- readr::read_delim
#
#  .read_block <- function(skip, nmax) {
#    read_delim(con, '\t',
#      col_names = c('CHR', 'SNPID', 'MORGANS', 'POS', 'A1', 'A2'),
#      col_types = readr::cols_only(SNPID = 'c'),
#      progress = FALSE,
#      skip = skip,
#      guess_max = 0,
#      n_max = nmax)
#  }
#
#  lst <- mapply(.read_block, skips, widths)
#
#
#}
