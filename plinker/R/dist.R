#' compute genotypic distances between samples
#'
#' @param type 	distance unit type.
#' 	cf <https://www.cog-genomics.org/plink/1.9/distance#read_dists>
#' @inheritParams bed_plink_cmd
#' @return a distance matrix
#'
#' @seealso bed_plink_distance
#' @export
bed_R_distance <- function(bo, type = c('ibs', 'allele-ct'))
{
  type <- match.arg(type)

  mat <- t(bed_genotypes(bo))
  nb <- bed_nb_samples(bo)

  dist <- matrix(NA, nrow = nb, ncol = nb)
  diag(dist) <- 1

  for (i in 1:(nb-1)) {
    vi <- mat[, i, drop = TRUE]
    for (j in (i+1):nb) {
      vj <- mat[, j, drop = TRUE]
      dist[j, i] <- dist[i, j] <- ibs(ibs_state(vi, vj))
    }
  }

  rownames(dist) <- colnames(dist) <- bed_sample_IDs(bo, custom = TRUE)

  dist
}

ibs <- function(state) {
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
