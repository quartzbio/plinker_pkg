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


#' run a plink command on a plinker_bed object
#'
#' @param ...		passed to \code{\link{plink_cmd}}
#' @inheritParams params
#' @inheritParams plink_cmd
#' @param stderr	cf \code{\link{system2}}
#' @export
bed_plink_cmd <- function(bo,
  args,
  snp_idx = bed_snp_idx(bo),
  sample_idx = bed_sample_idx(bo),
  quiet = FALSE,
  stdout = quiet,
  stderr = quiet,
  ...)
{
  ### add plink input files
  fns <- bo[c('bed', 'fam', 'bim')]
  fargs <- paste0('--', names(fns), ' ', fns)
  args <- c(fargs, args)

  ### add subsetting args if needed
  plink_snp_ids_fn <- 'input_snp_ids.txt'
  if (!is.null(snp_idx)) {
    bim_df <- bed_bim_df(bo, subset = FALSE)
    snp_ids <- bim_df$SNPID[snp_idx]

    writeLines(snp_ids, plink_snp_ids_fn)

    args <- c(paste('--extract', plink_snp_ids_fn), args)
  }

  plink_sample_ids_fn <- 'input_sample_ids.txt'
  if (!is.null(sample_idx)) {
    fam_df <- bed_fam_df(bo, subset = FALSE)
    sample_ids <- fam_df$FID[sample_idx]

    writeLines(sample_ids, plink_sample_ids_fn)

    args <- c(paste('--keep-fam', plink_sample_ids_fn), args)
  }

  invisible(plink_cmd(args, stdout = stdout, stderr = stderr, ...))
}


#' run plink --freq counts
#'
#' @param ...		passed to \code{\link{bed_plink_cmd}}
#' @inheritParams bed_plink_cmd
#' @return a data frame, cf \url{http://www.cog-genomics.org/plink/1.9/formats#frq_count}
#'
#' @export
bed_plink_freq_count <- function(bo, ...) {
  setup_temp_dir()
  bed_plink_cmd(bo, '--freq counts', ...)
  read_plink_freq_counts('plink.frq.counts')
}

#' get plink version
#'
#' @param ...		passed to \code{\link{plink_cmd}}
#' @export
plink_version <- function(...) {
  plink_cmd('--version', stdout = TRUE, ...)
}
