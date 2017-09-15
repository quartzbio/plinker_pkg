#' select a subsequence of samples in a plink dataset
#'
#' @inheritParams params
#' @return the subsetted bed dataset object
#' @family subset
#' @export
bed_subset_samples_by_idx <- function(bo, sample_idx) {
  if (length(sample_idx) == 0) stop('empty sample_idx')

  # check
  if (anyDuplicated(sample_idx) > 0) {
    stop('Error, you have duplicated sample indices')
  }

  if (any(is.na(sample_idx))) {
    stop('Error, you have missing (NA) sample indices')
  }

  # cleanup
  sample_idx <- as.integer(sample_idx)

  # check the range
  rg <- range(sample_idx)
  if (rg[1] < 1 || rg[2] > bed_nb_samples(bo))
    stop('bad sample_idx range')

  # apply the subset
  if (is.null(bo$sample_idx)) {
    # no subset yet
    bo$sample_idx <- sample_idx
  } else {
    bo$sample_idx <- bo$sample_idx[sample_idx]
  }

  bo
}


#' select a subset of samples in a plink dataset
#'
#' @inheritParams params
#' @return the subsetted bed dataset object
#' @family subset
#' @export
bed_subset_samples_by_IDs <- function(bo, sample_IDs) {
  idx <- bed_sample_IDs_to_idx(bo, sample_IDs)
  bed_subset_samples_by_idx(bo, idx)
}



#' cancel any subset of samples in a plink dataset
#'
#' @inheritParams params
#' @return the subsetted bed dataset object
#' @family subset
#' @export
bed_reset_subset_samples_by_idx <- function(bo) {
  bo$sample_idx <- NULL

  bo
}

#' convert sample ids to sample indices
#'
#' N.B: keep the ordering
#'
#' @inheritParams params
#' @return indices
#' @keywords internal
bed_sample_IDs_to_idx <- function(bo, sample_IDs, subset = TRUE,
  ignore_fid = bed_ignore_fid(bo))
{
  if (length(sample_IDs) == 0) stop('empty sample IDs')

  if (anyDuplicated(sample_IDs) > 0) {
    stop('Error, you have duplicated sample IDs')
  }

  all_ids <- bed_sample_IDs(bo, subset = subset, ignore_fid = ignore_fid)
  idx <- match(sample_IDs, all_ids)

  bads <- which(is.na(idx))
  if (length(bads) > 0) {
    bad_ids <- sample_IDs[bads]
    stop(sprintf('Error, bad sample IDs: "%s"', paste0(bad_ids, collapse = ',')))
  }


  idx
}
