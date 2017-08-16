#' select a subset of samples in a plink dataset
#'
#' @inheritParams params
#' @return the subsetted bed dataset object
#' @family subset
#' @export
bed_subset_samples_by_idx <- function(bo, sample_idx) {
  if (length(sample_idx) == 0) stop('empty sample_idx')

  # cleanup
  sample_idx <- sort(unique(as.integer(sample_idx)))

  # check the idx
  if (sample_idx[1] < 1 || sample_idx[length(sample_idx)] > bed_nb_samples(bo))
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
  current_ids <- bed_sample_IDs(bo, subset = TRUE)

  idx <- match(sample_IDs, current_ids)
  bads <- which(is.na(idx))
  if (length(bads) > 0) {
    bad_ids <- sample_IDs[bads]
    stop(sprintf('Error, bad sample IDs: "%s"', paste0(bad_ids, collapse = ',')))
  }

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
  all_ids <- bed_sample_IDs(bo, subset = subset, ignore_fid = ignore_fid)
  idx <- match(sample_IDs, all_ids)

  bads <- which(is.na(idx))
  if (length(bads) > 0) {
    bad_ids <- sample_IDs[bads]
    stop(sprintf('Error, bad sample IDs: "%s"', paste0(bad_ids, collapse = ',')))
  }


  idx
}
