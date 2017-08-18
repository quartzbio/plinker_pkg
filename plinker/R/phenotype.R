

#' make a phenotype vector to use with PLINK
#'
#'
#' Case/control phenotypes are expected to be encoded as:
#'   - 1=unaffected (control)
#'   - 2=affected (case)
#'   - -9=missing data
#'
#' cf <https://www.cog-genomics.org/plink/1.9/input#pheno>
#'
#' @inheritParams params
#' @param phenotype								the variable to consider
#' @param missing_phenotype 			integer value for replacing NA values in
#' 	the variable
#' @export
make_phenotype_from_vector <- function(phenotype, missing_phenotype = -9L)
{
  phenotype <- as.integer(phenotype)
  if (!is.integer(missing_phenotype) || length(missing_phenotype) != 1 ||
    is.na(missing_phenotype))
    stop('bad missing_phenotype: must integer of length 1')

  # check for occurrences of the missing_phenotype prior to substitute NA
  if (any(phenotype == missing_phenotype, na.rm = TRUE)) {
    stop('ERROR, the phenotype vector contain some "missing_phenotype"!')
  }

  phenotype[which(is.na(phenotype))] <- missing_phenotype

  phenotype
}


#' make a phenotype vector to use with PLINK from a data frame
#'
#' the data frame will be merged using [merge_df_with_fam()]
#'
#' @param	df							the data frame to extract the phenotype from
#' @param	phenotype_var		the name of the df column with the phenotype
#' @inheritDotParams merge_df_with_fam
#' @inheritParams params
#' @inheritParams make_phenotype_from_vector
#' @export
#' @seealso make_phenotype_from_vector
bed_phenotype_from_df <- function(bo,
  df,
  phenotype_var,
  missing_phenotype = -9L,
  subset = TRUE, ...)
{
  if (phenotype_var %in% df) stop('bad param phenotype_var')
  mdf <- merge_df_with_fam(bed_fam_df(bo, subset), df,
    ignore_fid = bed_ignore_fid(bo), ...)

  pheno <- make_phenotype_from_vector(mdf[[phenotype_var]], missing_phenotype)

  pheno
}
