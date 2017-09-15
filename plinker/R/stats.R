
#' compute linear regression using with covariates using R
#'
#' @inheritDotParams bed_plink_cmd -bo -args
#' @inheritParams bed_plink_lm
#' @export
#' @family plink

#' @seealso bed_phenotype_from_df
bed_R_lm <- function(bo,
  phenotype = bed_fam_df(bo)$PHENO,
  model = c('additive', 'dominant', 'recessive'),
  logistic = FALSE,
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

  if (logistic)
    phenotype <- as.factor(phenotype)

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

    fit <- if (logistic) {
      stats::glm(fo, data = df, model = TRUE, family = 'binomial')
    } else {
      stats::lm(fo, data = df, model = TRUE)
    }

    sm <- summary(fit)

    vars <- names(fit$coefficients[-1])
    res <- data.frame(
      bim[i, , drop = FALSE],
      TEST = vars,
      row.names = NULL,
      stringsAsFactors = FALSE)

    res$NMISS <- nrow(fit$model)

    # N.B: some vars may be absent from the summary
    if (logistic) {
      res$OR <- exp(fit$coefficients[vars])
      res$STAT <- sm$coefficients[, 'z value'][vars]
      res$P <- sm$coefficients[, 'Pr(>|z|)'][vars]
    } else {
      res$BETA <- fit$coefficients[vars]
      res$STAT <-sm$coefficients[, 't value'][vars]
      res$P <- sm$coefficients[, 'Pr(>|t|)'][vars]
    }

    res
  }

  dfs <- lapply(1:nb_snps, .process_snp)
  df <- do.call(rbind, dfs)

  df
}

#' compute fisher test using R
#'
#' currently only compute the genotypic fisher test
#'
#' @inheritParams params
#' @param phenotype		a binary phenotype vector as an integer vector
#' 	Case/control phenotypes are expected to be encoded as:
#'   - 1=unaffected (control)
#'   - 2=affected (case)
#' @export
#' @family stats
#' @seealso bed_phenotype_from_df
bed_R_fisher <- function(bo, phenotype = bed_fam_df(bo)$PHENO) {
  nb_snps <- bed_nb_snps(bo)

  bim <- bed_bim_df(bo)
  bim <- bim[, c('CHR', 'SNP', 'A1', 'A2')]

  .process_snp <- function(i) {
    genos <- bed_genotypes(bo, snp_idx = i)
    genos <- factor(genos, levels = 0:2)

    tt <- table(genos, phenotype)
    res <- stats::fisher.test(tt, conf.int = FALSE)

    df <- data.frame(
      bim[i, , drop = FALSE],
      TEST = 'GENO',
      row.names = NULL,
      AFF = paste0(tt[3:1,  2], collapse = '/'),
      UNAFF = paste0(tt[3:1,  1], collapse = '/'),
      P = res$p.value,
      stringsAsFactors = FALSE)

    df
  }

  dfs <- lapply(1:nb_snps, .process_snp)
  df <- do.call(rbind, dfs)

  df
}
