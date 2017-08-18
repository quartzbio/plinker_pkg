

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
#' @param variable								the variable to consider
#' @param missing_phenotype 			integer value for replacing NA values in
#' 	the variable
#' @export
#' @family phenotype
#' @md
bed_phenotype_from_vector <- function(bo,
  phenotype, missing_phenotype = -9L, subset = TRUE)
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
