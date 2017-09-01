
#' compute linear regression using with covariates using R
#'
#' @inheritDotParams bed_plink_cmd -bo -args
#' @inheritParams params
#' @inheritParams recode_genotypes
#' @param phenotype a quantitative phenotype/response vector
#' @param covars		the optional covars, as a data frame
#' @export
#' @family plink
#' @seealso make_phenotype_from_vector
#' @seealso bed_phenotype_from_df
bed_R_lm <- function(bo,
  phenotype = bed_fam_df(bo)$PHENO,
  model = c('additive', 'dominant', 'recessive'),
  covars = NULL
) {
  nb_snps <- bed_nb_snps(bo)
  model <- match.arg(model)

  var <- switch(model,
    additive = 'ADD',
    dominant = 'DOM',
    recessive = 'REC', stop('error'))

  bim <- bed_bim_df(bo)
  bim <- bim[, c('CHR', 'SNP', 'POS', 'A1')]
  names(bim) <- c('CHR', 'SNP', 'BP', 'A1')

  .process_snp <- function(i) {
    genos <- bed_genotypes(bo, snp_idx = i)

    X <- as.numeric(recode_genotypes(genos, model))
    df <- data.frame(phenotype, X, stringsAsFactors = FALSE)
    names(df) <- c('Y', var)
    fo <- sprintf('Y ~ %s', var)
    if (!is.null(covars)) {
      cdf <- covars[, -(1:2), drop = FALSE]
      fo <- paste0(fo, ' + ', paste0(names(cdf), collapse = ' + '))
      df <- cbind(df, cdf)
    }

    fit <- stats::lm(fo, data = df, model = TRUE)
    sm <- summary(fit)

    vars <- names(fit$coefficients[-1])
    res <- data.frame(
      bim[i, , drop = FALSE],
      TEST = vars,
      row.names = NULL,
      stringsAsFactors = FALSE)

    res$NMISS <- nrow(fit$model)
    # N.B: some vars may be absent from the summary
    res$BETA <- fit$coefficients[vars]
    res$STAT <-sm$coefficients[, 't value'][vars]
    res$P <- sm$coefficients[, 'Pr(>|t|)'][vars]

    res
  }

  dfs <- lapply(1:nb_snps, .process_snp)
  df <- do.call(rbind, dfs)

  df
}

