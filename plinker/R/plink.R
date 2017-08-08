#' run a plink command using system2
#'
#' @param ...		passed to \code{\link{system2}}
#' @inheritParams base::system2
#' @keywords internal
plink_cmd <- function(args, command = Sys.which('plink'),
  stdout = "", stderr = "", ...)
{
  res <- system2(command, args, stdout = stdout, stderr = stderr, ...)

  status <- res
  if (isTRUE(stdout) || isTRUE(stderr)) {
    status <- attr(res, 'status')
    if (is.null(status)) status <- 0
  }

  if (status != 0) {
    cmd <- paste(command, args, collapse = ' ')
    msg <- sprintf('ERROR running plink "%s", exit status=%i', cmd, status)
    stop(msg)
  }

  res
}

#' get plink version
#'
#' @param ...		passed to \code{\link{plink_cmd}}
#' @export
plink_version <- function(...) {
  plink_cmd('--version', stdout = TRUE, ...)
}
