
#' group a set of integers by contiguous blocks
#'
#' @param ints		the integers as a sorted vector
#' @return a 2 column
#' @keywords internal
split_sorted_ints_by_blocks <- function(ints) {
  nb <- length(ints)
  if (nb == 0) return(NULL)

  # try to be smart by testing obvious use cases first
  if (ints[nb] == ints[1] + nb - 1) {
    return(cbind(ints[1], nb, deparse.level = 0))
  }

  starts <- c(1L, which(diff(ints) != 1L) + 1L)
  lens <- diff(c(starts, nb + 1L))
  blks <- cbind(ints[starts], lens, deparse.level = 0)

  blks
}
