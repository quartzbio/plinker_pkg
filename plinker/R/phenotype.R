

check_missing_phenotype <- function(missing_phenotype) {
  if (!is.integer(missing_phenotype) || length(missing_phenotype) != 1 ||
    is.na(missing_phenotype) || missing_phenotype >= 0)
    stop('bad missing_phenotype: must be negative integer of length 1')
}

check_case_control_phenotype <- function(phenotype, nb_samples) {
  if (length(phenotype) != nb_samples)
    stop('bad phenotype length')
  if (!is.integer(phenotype)) stop('phenotype must be an integer vector!')

  values <- stats::na.omit(unique(phenotype))
  if (length(values) > 2) stop('phenotype must be binary')
}


check_quantitative_phenotype <- function(phenotype, nb_samples) {
  if (length(phenotype) != nb_samples)
    stop('bad phenotype length')
  if (!is.numeric(phenotype)) stop('phenotype must be a numeric vector!')
}


#' make a phenotype vector to use with PLINK from a data frame
#'
#' the data frame will be merged using [merge_df_with_fam()]
#'
#' @param	df							the data frame to extract the phenotype from
#' @param	phenotype_var		the name of the df column with the phenotype
#' @inheritDotParams merge_df_with_fam
#' @inheritParams params
#' @export
bed_phenotype_from_df <- function(bo,
  df,
  phenotype_var,
  subset = TRUE, ...)
{
  if (phenotype_var %in% df) stop('bad param phenotype_var')
  mdf <- merge_df_with_fam(bed_fam_df(bo, subset), df,
    ignore_fid = bed_ignore_fid(bo), ...)

  mdf[[phenotype_var]]
}
