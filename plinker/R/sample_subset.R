#' select a subset of samples in a plink dataset
#'
#' @inheritParams params
#' @return the subset bed dataset object
#' @family accessors
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


#' cancel any subset of samples in a plink dataset
#'
#' @inheritParams params
#' @return the unsubset bed dataset object
#' @family accessors
#' @export
bed_reset_subset_samples_by_idx <- function(bo) {
  bo$sample_idx <- NULL

  bo
}

