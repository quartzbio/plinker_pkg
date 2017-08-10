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
#' @param allow_no_sex		if FALSE, samples with ambiguous sex have their
#' 		 phenotypes set to missing, cf the \code{--allow-no-sex} option in
#' 			\url{http://www.cog-genomics.org/plink/1.9/filter#maf}
#' @param nonfounders		by default in plink, nonfounders are not counted
#' 		by --freq{x} or --maf/--max-maf/--hwe. setting this to TRUE disable
#' 		this behaviour by using the \code{--nonfounders} option, cf
#' 		\url{http://www.cog-genomics.org/plink/1.9/filter#nonfounders}
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
  allow_no_sex = TRUE,
  nonfounders = TRUE,
  quiet = FALSE,
  stdout = if (quiet) FALSE else "",
  stderr = if (quiet) FALSE else "",
  ...)
{
  ### add plink input files
  fns <- bo[c('bed', 'fam', 'bim')]
  fargs <- paste0('--', names(fns), ' ', fns)
  args <- c(fargs, args)

  ## add optional flags
  if (allow_no_sex) args <- c(args, '--allow-no-sex')
  if (nonfounders) args <- c(args, '--nonfounders')
  ### add subsetting args if needed
  plink_snp_ids_fn <- 'input_snp_ids.txt'
  if (!is.null(snp_idx)) {
    bim_df <- bed_bim_df(bo, subset = FALSE)
    snp_ids <- bim_df$SNP[snp_idx]

    writeLines(snp_ids, plink_snp_ids_fn)

    args <- c(paste('--extract', plink_snp_ids_fn), args)
  }

  plink_sample_ids_fn <- 'input_sample_ids.txt'
  if (!is.null(sample_idx)) {
    fam_df <- bed_fam_df(bo, subset = FALSE)

    readr::write_tsv(fam_df[sample_idx, 1:2, drop = FALSE],
      plink_sample_ids_fn, col_names = FALSE)

    args <- c(paste('--keep', plink_sample_ids_fn), args)
  }

  invisible(plink_cmd(args, stdout = stdout, stderr = stderr, ...))
}


#' compute allele frequencies using plink --freq counts
#'
#' @param ...		passed to \code{\link{bed_plink_cmd}}
#' @inheritParams bed_plink_cmd
#' @return a data frame, cf \url{http://www.cog-genomics.org/plink/1.9/formats#frq_count}
#'
#' @seealso bed_plink_cmd
#' @export
bed_plink_freq_count <- function(bo, ...)
 {
  setup_temp_dir()

  args <- '--freq counts'
  bed_plink_cmd(bo, args, ...)
  read_plink_freq_counts('plink.frq.counts')
}

#' compute genotype frequencies using plink --freqx
#'
#' @param ...		passed to \code{\link{bed_plink_cmd}}
#' @inheritParams bed_plink_cmd
#' @return a data frame, cf \url{http://www.cog-genomics.org/plink/1.9/formats#frqx}
#'
#' @seealso bed_plink_cmd
#' @export
bed_plink_freqx <- function(bo, ...)
{
  setup_temp_dir()

  args <- '--freqx'
  bed_plink_cmd(bo, args, ...)
  read_plink_output('plink.frqx')
}

#' compute snp and sample missing rates using plink --missing
#'
#' @param ...		passed to \code{\link{bed_plink_cmd}}
#' @inheritParams bed_plink_cmd
#' @return a list of 2 data frame,
#' 	cf \url{http://www.cog-genomics.org/plink/1.9/formats#imiss} and
#' 	\url{http://www.cog-genomics.org/plink/1.9/formats#lmiss}
#'
#' @seealso bed_plink_cmd
#' @export
bed_plink_missing <- function(bo, ...)
{
  setup_temp_dir()

  args <- '--missing'
  bed_plink_cmd(bo, args, ...)
  imiss <- read_plink_output('plink.imiss')
  lmiss <- read_plink_output('plink.lmiss')

  list(sample = imiss, snp = lmiss)
}


#' get plink version
#'
#' @param ...		passed to \code{\link{plink_cmd}}
#' @export
plink_version <- function(...) {
  plink_cmd('--version', stdout = TRUE, ...)
}
