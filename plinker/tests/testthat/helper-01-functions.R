
fork_it <- function(..., silent = FALSE) {
  ### N.B: we use explicit parallel:: to be able to use this
  ### function when bootstrapping qbutils (i.e. no doc or outdated doc)
  res <- parallel::mccollect(parallel::mcparallel(silent = silent, ...))[[1]]

  if (is_error(res)) stop(get_error(res))

  invisible(res)
}


is_error <- function(x) {
  inherits(x, 'try-error')
}


get_error <- function(x) {
  if (!is_error(x)) stop('x must be an error')
  attr(x, 'condition')
}


get_error_msg <- function(x) {
  conditionMessage(get_error(x))
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

