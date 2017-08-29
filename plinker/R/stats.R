
#' compute linear regression using with covariates using R
#'
#' @inheritDotParams bed_plink_cmd -bo -args
#' @inheritParams params
#' @inheritParams recode_genotypes
#' @param phenotype a quantitative phenotype/response vector
#' @export
#' @family plink
#' @seealso make_phenotype_from_vector
#' @seealso bed_phenotype_from_df
bed_R_lm <- function(bo,
  phenotype = bed_fam_df(bo)$PHENO,
  model = c('additive', 'dominant', 'recessive'))
{
  browser()

  nb_snps <- bed_nb_snps(bo)

  .process_snp <- function(i) {
    genos <- bed_genotypes(bo, snp_idx = i)
    X <- recode_genotypes(genos, model)
    df <- data.frame(X = X, Y = phenotype, stringsAsFactors = FALSE)
    fit <- stats::lm(Y ~ X, data = df)
  }


}

