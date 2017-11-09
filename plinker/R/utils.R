
#' group a set of integers by contiguous blocks
#'
#' @param ints		the integers as a sorted vector
#' @return a 2 column
#' @keywords internal
split_sorted_ints_by_blocks <- function(ints) {
  ints <- as.integer(ints)
  nb <- length(ints)
  if (nb == 0) return(NULL)

  # try to be smart by testing obvious use cases first
  if (ints[nb] == ints[1] + nb - 1) {
    return(cbind(ints[1], nb, deparse.level = 0))
  }

  starts <- c(1L, which(diff(ints) != 1L) + 1L)
  lens <- diff(c(starts, nb + 1L))
  blks <- cbind(ints[starts], lens, deparse.level = 0)

  blks
}



# setup a temp dir
setup_temp_dir <- function() {
  caller <- as.character(sys.call(-1))[1]

  dir <- tempfile(caller)
  dir.create(dir, recursive = TRUE)

  old <- setwd(dir)
  # just to avoid warnings when tested via qbdev
  writeLines('', '.qbdev')

  cleanup <- bquote({
      setwd(.(old))
      unlink(.(dir), recursive = TRUE)
    })
  env <- parent.frame()

  do.call(add_on_exit, list(cleanup, env))

  invisible(dir)
}

add_on_exit <- function(expr, where = parent.frame()) {
  do.call("on.exit", list(substitute(expr), add = TRUE), envir = where)
}

