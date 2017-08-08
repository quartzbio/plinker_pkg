#' select a subset of SNPs in a plink dataset
#'
#' @inheritParams params
#' @return the subset bed dataset object
#' @family accessors
#' @export
bed_subset_snps_by_idx <- function(bo, snp_idx) {
  if (length(snp_idx) == 0) stop('empty snp_idx')

  # cleanup
  snp_idx <- sort(unique(as.integer(snp_idx)))

  # check the idx
  if (snp_idx[1] < 1 || snp_idx[length(snp_idx)] > bed_nb_snps(bo))
    stop('bad snp_idx range')

  # apply the subset
  if (is.null(bo$snp_idx)) {
    # no subset yet
    bo$snp_idx <- snp_idx
  } else {
    bo$snp_idx <- bo$snp_idx[snp_idx]
  }

  bo
}


#' cancel any subset of SNPs in a plink dataset
#'
#' @inheritParams params
#' @return the unsubset bed dataset object
#' @family accessors
#' @export
bed_reset_subset_snps_by_idx <- function(bo) {
  bo$snp_idx <- NULL

  bo
}
