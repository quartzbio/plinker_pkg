#' compute genotypic distances between samples
#'
#' @param type 	distance unit type. Currently only ibs
#' 	cf <https://www.cog-genomics.org/plink/1.9/distance#read_dists>
#' @inheritSection ibs  missingness
#' @inheritParams bed_plink_cmd
#' @return a distance matrix
#'
#' @seealso bed_plink_distance
#' @export
bed_R_distance <- function(bo, type = c('ibs'))
{
  type <- match.arg(type)

  dist <- compute_ibs_matrix(bed_genotypes(bo))
  rownames(dist) <- colnames(dist) <- bed_sample_IDs(bo, custom = TRUE)

  dist
}

compute_ibs_matrix <- function(genotypes) {
  mat <- t(genotypes)
  nb <- ncol(mat)

  dist <- matrix(NA, nrow = nb, ncol = nb)
  diag(dist) <- 1

  for (i in 1:(nb-1)) {
    vi <- mat[, i, drop = TRUE]
    for (j in (i+1):nb) {
      vj <- mat[, j, drop = TRUE]
      dist[j, i] <- dist[i, j] <- ibs(vi, vj)
    }
  }

  dist
}


#' compute the IBS similarity between two vector of genotypes
#'
#' Beware: the lengths of genotypes vectors are not checked
#'
#' @section missingness:
#' Currently, it uses what PLINK calls _flat-missing_ missingness correction,
#' which is the distance when only the non-missing pairs of genotypes
#' are considered
#'
#' @param genos1	first vector of genotypes (NA, 0:2)
#' @param genos2	second vector of genotypes (NA, 0:2)
#' @param state		the IBS state vector on which to compute the IBS
#' @return the IBS distance as a numeric scalar, or NA if no genotype to compare
#'
#' @export
ibs <- function(genos1, genos2, state = ibs_state(genos1, genos2)) {
  nas <- sum(is.na(state))
  nb <- length(state) - nas
  if (nb >= 0)
      sum(state, na.rm = TRUE) / (nb + nb)
  else
    NA_real_
}

ibs_state <- function(genos1, genos2) {
  2L - abs(genos1 - genos2)
}
