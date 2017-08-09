
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

setup_temp_dir <- plinker:::setup_temp_dir





