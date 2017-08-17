
#' convert a genotype matrix of raw genotypes to strings
#'
#' The genotypes must be encoded as [bed_genotypes]
#'
#' @param genos	a genotype matrix, as returned by [bed_genotypes]
#' @param allele1		the first allele, as a character vector of length the
#' 		nb of columns of genos
#' @param allele2		the second allele, as a character vector of length the
#' 		nb of columns of genos
#' @param sep				the separator for the genotype strings, as a character scalar
#' @param	sort			whether to sort the (heterozygote) strings
#' @return a character matrix of the same dimension as genos
#' @export
#' @md
convert_genotypes_to_string <- function(genos,
  allele1,
  allele2,
  sep = '/',
  sort = TRUE)
{

  ## check args
  if (!is.matrix(genos) || !is.integer(genos))
    stop('bad param genos: must be a matrix of integers')
  nsnps <- ncol(genos)

  if (!is.character(allele1) || length(allele1) != nsnps)
    stop('bad param allele1')
  if (!is.character(allele2) || length(allele2) != nsnps)
    stop('bad param allele2')
  if (!is.character(sep) || length(sep) != 1)
    stop('bad param sep')

  df <- make_genotype_converters(allele1, allele2, sep, sort)
  convs <- t(df[, 4:6])
  colnames(convs) <- df$type
  rownames(convs) <- NULL

  types <- paste0(allele1, sep, allele2)

  .convert_col <- function(i) {
    type <- types[i]
    conv <- convs[, type, drop = TRUE]
    strs <- conv[genos[, i, drop = TRUE] + 1L] # map to vector indices
  }

  mat <- sapply(seq_len(nsnps), .convert_col)
  dimnames(mat) <- dimnames(genos)


  mat
}


make_genotype_converters <- function(allele1, allele2, sep = '/', sort = TRUE) {
  types <- unique(paste0(allele1, sep, allele2))

  df <- data.frame(A1 = allele1, A2 = allele2,
    type = paste0(allele1, sep, allele2), stringsAsFactors = FALSE)

  df <- unique(df)

  a1 <- df$A1
  a2 <- df$A2
  a11 <- a1
  a22 <- a2
  if (sort) {
    a11 <- pmin(a1, a2)
    a22 <- pmax(a1, a2)
  }

  conv <- cbind(
      A2A2 = paste0(a2, sep, a2),
      A1A2 = paste0(a11, sep, a22),
      A1A1 = paste0(a1, sep, a1))

  res <- cbind(df, conv, stringsAsFactors = FALSE)

  res
}
