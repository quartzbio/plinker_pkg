#' select a subset of SNPs in a plink dataset
#'
#' N.B: the order of SNPs will not be preserved.
#'
#' @inheritParams params
#' @return the subsetted bed dataset object
#' @family subset
#' @export
bed_subset_snps_by_idx <- function(bo, snp_idx) {
  if (length(snp_idx) == 0) stop('empty snp_idx')

  if (anyDuplicated(snp_idx) > 0)
    stop('duplicated snp indices')

  if (any(is.na(snp_idx))) {
    stop('Error, you have missing (NA) snp indices')
  }

  snp_idx <- as.integer(snp_idx)

  # check the idx
  rg <- range(snp_idx)
  if (rg[1] < 1 || rg[2] > bed_nb_snps(bo))
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

#' select a subset of SNPs in a plink dataset
#'
#' @inheritParams params
#' @return the subsetted bed dataset object
#' @family subset
#' @export
bed_subset_snps_by_IDs <- function(bo, snp_IDs) {
  bim_df <- bed_bim_df(bo, subset = TRUE)

  idx <- match(snp_IDs, bim_df$SNP)
  bads <- which(is.na(idx))
  if (length(bads) > 0) {
    bad_ids <- snp_IDs[bads]
    stop(sprintf('Error, bad snp_IDs: "%s"', paste0(bad_ids, collapse = ',')))
  }

  bed_subset_snps_by_idx(bo, idx)
}


#' cancel any subset of SNPs in a plink dataset
#'
#' @inheritParams params
#' @return the unsubsetted bed dataset object
#' @family subset
#' @export
bed_reset_subset_snps_by_idx <- function(bo) {
  bo$snp_idx <- NULL

  bo
}
