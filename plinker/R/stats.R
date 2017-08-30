
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
  nb_snps <- bed_nb_snps(bo)

  model <- match.arg(model)

  test <- switch(model,
    additive = 'ADD',
    dominant = 'DOM',
    recessive = 'REC', stop('error'))

  .process_snp <- function(i) {
    genos <- bed_genotypes(bo, snp_idx = i)

    X <- recode_genotypes(genos, model)
    df <- data.frame(X = X, Y = phenotype, stringsAsFactors = FALSE)
    fit <- stats::lm(Y ~ X, data = df, model = FALSE)
    sm <- summary(fit)

    data.frame(
      NMISS = sum(!is.na(X)),
      BETA = fit$coefficients[2],
      STAT = sm$coefficients['X', 't value'],
      P = sm$coefficients['X', 'Pr(>|t|)'],
      stringsAsFactors = FALSE)
  }

  dfs <- lapply(1:nb_snps, .process_snp)
  df <- do.call(rbind, dfs)

  bim <- bed_bim_df(bo)
  bim <- bim[, c('CHR', 'SNP', 'POS', 'A1')]
  names(bim) <- c('CHR', 'SNP', 'BP', 'A1')

  cbind(bim, df)
}

