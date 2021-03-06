
#' init bedmatrix if needed
#'
#' @inheritParams params
#' @return an initialized BEDMatrix object
#' @keywords internal
bed_init_bedmatrix <- function(bo) {
  # is it already defined ?
#  if (!is.null(bo$cache$bedmatrix)) return(invisible())

  mat <- BEDMatrix::BEDMatrix(bo$bed,
    n = bed_nb_samples(bo, subset = FALSE),
    p = bed_nb_snps(bo ,subset = FALSE))

  # set dimnames
  dimnames(mat) <- list(bed_sample_IDs(bo, subset = FALSE),
    bed_snp_IDs(bo, subset = FALSE))

#  bo$cache$bedmatrix <- mat

  mat
}

#' get genotypes as an integer matrix
#'
#' The genotypes are coded as in plink **.raw** files:
#'   - 0 indicates A2A2
#'   - 1 indicates A1A2
#'   - 2 indicates A1A1
#'   - NA indicates missing
#'
#' @inheritParams params
#' @param subset	whether to consider subset info from the bo object.
#' @param lexicographic_allele_order
#' 								whether to use lexicographic order for alleles, and thus
#' 								to recode genotypes accordingly
#' @return the genotypes as an integer matrix of samples X snps,
#' @export
#' @md
bed_genotypes <- function(bo, subset = TRUE, lexicographic_allele_order = FALSE) {
  bmat <- bo$cache$bedmatrix
  if (is.null(bmat))
    bmat <- bo$cache$bedmatrix <- bed_init_bedmatrix(bo)

  snp_idx <- bed_snp_idx(bo)
  sample_idx <- bed_sample_idx(bo)


  genos <- NULL
  if (!subset || (is.null(snp_idx) && is.null(sample_idx))) {
    # N.B: the empty [] calls BedMatrix method
    genos <- bmat[,, drop = FALSE]
  } else {
    genos <- if (is.null(snp_idx)) {
      bmat[sample_idx, , drop = FALSE]
    } else if (is.null(sample_idx)) {
      bmat[, snp_idx, drop = FALSE]
    } else {
      bmat[sample_idx, snp_idx, drop = FALSE]
    }
  }

  if (lexicographic_allele_order) {
    a2h <- bed_allele_higher(bo, subset = subset)
    a2 <- bed_allele2(bo, subset = subset)
    inv <- which(a2h != a2)

    genos[, inv] <- 2L - genos[, inv]
  }

  genos
}


