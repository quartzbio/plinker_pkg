
#' build a long data frame with useful information about the genotypes
#'
#' **Warning**: very memory intensive, only use on small datasets/subsets
#'
#' @keywords internal
bed_convert_genotypes_to_data_frame <- function(bo, subset = TRUE) {
  genos <- bed_genotypes(bo, subset)
  genos_str <- bed_genotypes_as_strings(bo, subset)

  bim <- bed_bim_df(bo, subset)
  fam <- bed_fam_df(bo, subset)

  sample_ids <- rownames(genos)

  .snp_to_df <- function(i) {
    gi <- genos[, i]
    data.frame(
      bim[i, ],
      SAMPLE_ID = sample_ids,
      GENO_INT = gi,
      GENO_STR = genos_str[, i],
      DOMINANT = recode_genotypes(gi, 'dominant'),
      RECESSIVE = recode_genotypes(gi, 'recessive'),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  dfs <- lapply(1:ncol(genos), .snp_to_df)

  do.call(rbind.data.frame, dfs)
}
