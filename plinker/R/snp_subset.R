#' select a subset of SNPs in a plink dataset
#'
#' N.B: the order of SNPs will not be preserved.
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

#' select a subset of SNPs in a plink dataset
#'
#' @inheritParams params
#' @return the subset bed dataset object
#' @family accessors
#' @export
bed_subset_snps_by_IDs <- function(bo, snp_IDs) {
  bim_df <- bed_bim_df(bo, subset = TRUE)

  snp_IDs <- unique(snp_IDs)
  idx <- match(snp_IDs, bim_df$SNPID)
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
#' @return the unsubset bed dataset object
#' @family accessors
#' @export
bed_reset_subset_snps_by_idx <- function(bo) {
  bo$snp_idx <- NULL

  bo
}
