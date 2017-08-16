
#' read and parse a plink .fam file
#'
#' cf \url{http://www.cog-genomics.org/plink/1.9/formats#fam}
#'
#' @param path		the .fam file path
#' @return a data frame
#' @export
read_fam <- function(path) {
  as.data.frame(readr::read_delim(path, ' ',
    col_names = c('FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'SEX', 'PHENO'),
    col_types = 'ccccii',
    progress = FALSE
   ))
}


