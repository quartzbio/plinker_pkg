

#' run a plink command using system2
#'
#' @inheritParams params
#' @param ...		passed to \code{\link{system2}}
#' @inheritParams base::system2
#' @keywords internal
plink_cmd <- function(args, command = Sys.which('plink'),
  quiet = FALSE,
  stdout = if (quiet) FALSE else "",
  stderr = if (quiet) FALSE else "",
  ...)
{
  res <- system2(command, args, stdout = stdout, stderr = stderr, ...)

  status <- res
  if (isTRUE(stdout) || isTRUE(stderr)) {
    status <- attr(res, 'status')
    if (is.null(status)) status <- 0
  }

  if (status != 0) {
    args <- paste(args, collapse = ' ')
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


#' convert a PLINK transposed tped dataset (.tped, .tfam) to a bed dataset
#'
#' PLINK will decide which alleles are A1 and A2 based on allele frequency.
#' A1 will be the minor allele, and on tie use the lexicographic order.
#'
#' @param input_prefix		the prefix for the .tped/.tfam input files
#' @param output_prefix   the prefix for the output (.bed, .bim, .fam) files
#' @param not_chr 				chromosomes to exclude
#' @inheritDotParams plink_cmd -args
#' @export
#' @seealso plink_bfile_to_tfile
plink_tfile_to_bfile <- function(input_prefix, output_prefix, not_chr = '0', ...) {

  args <- c(
    sprintf('--tfile %s', input_prefix),
    sprintf('--out %s', output_prefix),
    sprintf('--not-chr %s', not_chr)
    )

  plink_cmd(args, ...)
}

#' save a new PLINK binary dataset with alleles ordered lexicographically
#'
#' @param input_prefix		the prefix for the input bed dataset
#' @param output_prefix   the prefix for the output bed dataset
#' @inheritDotParams plink_cmd -args
#' @export
save_plink_with_lexicographic_alleles_order <- function(input_prefix,
  output_prefix, ...) {

  bo <- bed_open(input_prefix)
  bim <- bed_bim_df(bo)
  bim$a2max <- with(bim, pmax(A1, A2))

  alleles_file <- tempfile()
  save_bim(bim[, c('SNP', 'a2max')], alleles_file)
  on.exit(unlink(alleles_file), add = TRUE)

  args <- c(
    sprintf('--bfile %s', input_prefix),
    sprintf('--a2-allele  %s', alleles_file),
    '--make-bed',
    sprintf('--out %s', output_prefix)
  )

  plink_cmd(args, ...)
}


#' convert a  a bed dataset dataset to a PLINK transposed tped dataset
#'
#'
#' @param input_prefix		the prefix for the output (.bed, .bim, .fam) files
#' @param output_prefix		the prefix for the output .tped/.tfam input files
#' @inheritDotParams plink_cmd -args
#' @export
#' @seealso plink_tfile_to_bfile
plink_bfile_to_tfile <- function(input_prefix, output_prefix, ...) {
  args <- c(
    sprintf('--bfile %s', input_prefix),
    sprintf('--recode transpose --out %s', output_prefix),
    '--allow-no-sex',
    '--nonfounders',
    '--keep-allele-order'
#      ,sprintf('--maf %s', maf)
  )
  plink_cmd(args, ...)
}

#' write a PLINK alleles file, as used by --update-alleles or --recoded-alleles
#'
#' @param alleles_df		 a data frame, whose 3 columns are the variant ID, the allele 1
#' 	and the allele 2
#' @param path	the output path
#' @family write
#' @export
save_plink_alleles <- function(alleles_df, path) {
  if (!is.data.frame(alleles_df) && ncol(alleles_df) >= 3)
    stop('bad arg alleles_df')

  utils::write.table(alleles_df, path, sep = ' ',
     col.names = FALSE, row.names = FALSE, quote = FALSE)
}

