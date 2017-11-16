
#' filter snps from a bed by maf
#'
#' @inheritParams params
#' @param max		the maximum snp missing rate
#' @param ...		passed to \code{\link{bed_plink_cmd}}
#'
#' @return the subsetted bed dataset object
#' @family accessors
#' @export
bed_filter_snps_by_missing_rate <- function(bo, max, ...) {
  df <- bed_plink_missing(bo, ...)$snp
  goods <- which(df$F_MISS <= max)

  bed_subset_snps_by_idx(bo, goods)
}
