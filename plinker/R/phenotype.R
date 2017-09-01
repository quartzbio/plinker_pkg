

#' make a phenotype vector to use with PLINK
#'
#'
#' Case/control phenotypes are expected to be encoded as:
#'   - 1=unaffected (control)
#'   - 2=affected (case)
#'
#' Missing phenotype is encoded by -9 by default.
#'
#'
#' cf <https://www.cog-genomics.org/plink/1.9/input#pheno>
#'
#' @inheritParams params
#' @param phenotype								the variable to consider
#' @param missing_phenotype 			integer value for replacing NA values in
#' 																the variable
#' @export
make_phenotype_from_vector <- function(phenotype, missing_phenotype = -9L)
{
  nas_ind <- which(is.na(phenotype))
  if (length(nas_ind) == 0) return(phenotype)

  check_missing_phenotype(missing_phenotype)

  # check for occurrences of the missing_phenotype prior to substitute NA
  if (any(phenotype == missing_phenotype, na.rm = TRUE)) {
    stop('ERROR, the phenotype vector contain some "missing_phenotype"!')
  }

  phenotype[nas_ind] <- missing_phenotype

  phenotype
}

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
#' @seealso make_phenotype_from_vector
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
