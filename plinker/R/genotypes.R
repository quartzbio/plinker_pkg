#' get genotypes as a matrix of strings
#'
#'
#' @inheritParams params
#' @inheritParams convert_genotypes_to_string
#' @param subset	whether to consider subset info from the bo object.
#' 	if FALSE, no other params must be given
#' @return the genotypes as an integer matrix of samples X snps,
#' @export
bed_genotypes_as_strings <- function(bo, subset = TRUE, sep = '/', sort = TRUE)
{
  allele1 <- bed_allele1(bo, subset)
  allele2 <- bed_allele2(bo, subset)
  mat <- bed_genotypes(bo, subset = subset)
  strs <- convert_genotypes_to_string(mat, allele1, allele2, sep = sep,
    sort = sort)

  strs
}


#' convert a genotype matrix of raw genotypes to strings
#'
#' The genotypes must be encoded as [bed_genotypes()]
#'
#' @param genos	a genotype matrix, as returned by [bed_genotypes()]
#' @param allele1		the first allele, as a character vector of length the
#' 		nb of columns of genos
#' @param allele2		the second allele, as a character vector of length the
#' 		nb of columns of genos
#' @param sep				the separator for the genotype strings, as a character scalar
#' @param	sort			whether to sort the (heterozygote) strings
#' @return a character matrix of the same dimension as genos
#' @export
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
  # special case needed when mat is a scalar
  dim(mat) <- dim(genos)
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

extract_genotypes_from_ped <- function(ped_df, sep = '/') {
  mat <- as.matrix(ped_df[, -(1:6)])
  # paste all pairs of columns

  .paste_pair <- function(i) {
    paste0(mat[, i], '/', mat[, i + 1])
  }
  strs2 <- sapply(seq.int(start = 1, length.out = ncol(mat)/2, by = 2),
    .paste_pair)
  strs2[strs2 == '0/0'] <- NA

  strs2
}


#' recode numeric genotype according to a genetic model
#'
#' N.B: in PLINK, the A1 allele is usually the minor allele, but that
#' is often false. Moreover if you subset your dataset, it may invert the
#' minor and major allele.
#'
#' Here all the models assume that A1 is the risk/trait allele. You can use
#' the invert param to use A2 instead of A1.
#'
#' @param genos 		a numeric vector of genotypes (for one SNP). The values are
#'	the same as [bed_genotypes()], i.e. A2A2=0, A1A2=1, A1A1=2,NA=missing
#'
#' @param model 	the genetic model either:
#' 		- additive: (A2A2=0, A1A2=1, A1A1=2), i.e the the number of A1 alleles
#' 				i.e. no change
#'    - dominant: (A2A2=0, A1A2=1, A1A1=1), i.e presence of A1
#'    - recessive: (A2A2=0, A1A2=0, A1A1=1), i.e. presence of A1A1
#'
#' @param invert	whether to consider A2 as the risk/trait allele
#'
#' @return  numeric vector of recoded genotypes
#'
#' @family genotypes
#' @keywords internal
recode_genotypes <- function(genos,
  model = c('additive', 'dominant', 'recessive'), invert = FALSE)
{
  model <- match.arg(model)

  if (invert)
    genos <- 2L - genos

  switch(model,
    additive = genos,
    dominant = pmin(genos, 1L),
    recessive = pmax(genos - 1L, 0L),
    stop("not possible")
  )
}





