

#' create a bed object
#'
#' @param bed		the .bed filename
#' @param bim		the .bim
#' @keywords internal
new_bed <- function(bed, bim, fam) {

  # user friendly check
  files <- list(bed = bed, bim = bim, fam = fam)
  for (i in seq_along(files)) {
    if (!file.exists(files[[i]])) {
      msg <- sprintf('Error, %s file "%s" does not exist.', names(files)[i],
        files[[i]])
      stop(msg)
    }
  }

  obj <- normalizePath(unlist(files))
  names(obj) <- names(files)
  obj <- as.list(obj)

  class(obj) <- 'plinker_bed'

  obj
}



unprefix_bed <- function(prefix) {
  types <- c('bed', 'bim', 'fam')
  fns <- paste0(prefix, '.', types)
  names(fns) <- types

  fns
}
