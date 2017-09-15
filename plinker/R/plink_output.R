#' read a plink ped file
#'
#' @param path	the path of the ped file
#' @inheritParams params
#' @return the parsed ped as a data frame
#' @family read
#' @export
read_plink_ped <- function(path) {
  df <- read_plink_output(path)
  nsnps <- (ncol(df) - 6) / 2

  snp_cols <- rep(paste0('SNP', seq_len(nsnps)), each = 2)
  snp_cols <- paste0(snp_cols, c('A1', 'A2'))

  names(df) <- c('FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'SEX', 'PHENO', snp_cols)

  df
}

read_plink_output <- function(path, ...) {
  if (!file.exists(path))
    stop('plink output file does not exist:', path)


  df <- data.table::fread(path,
    verbose = FALSE,
    data.table = FALSE,
    showProgress = FALSE,
    ...
  )

  if (!is.null(df$CHR)) df$CHR <- as.character(df$CHR)

  df
}

read_plink_freq <- function(path) {
  read_plink_output(path)
  #read_plink_output(path, col_types = 'ccccni')
}

read_plink_freq_counts <- function(path) {
  read_plink_output(path,
    colClasses = c('character','character', 'character', 'character',
      'integer', 'integer', 'integer'))
#  read_plink_output(path, col_types = 'cccciii')
}

#' reorder if needed a plink output
#'
#' useful to apply the ordering of a bed_plink object
#'
#' @param df		a PLINK output as a data frame, that has a SNP column
#' @inheritParams params
#' @return the reordered output as a data frame
#' @keywords internal
reorder_plink_output <- function(df, snp_IDs) {
  if (!'SNP' %in% names(df)) stop('no SNP column')

  if (length(snp_IDs) == 0) return(df)

  if (!all(df$SNP %in% snp_IDs))  stop('bad SNP ids')

  snp_IDs <- intersect(snp_IDs, df$SNP)

  tt <- table(df$SNP)
  tt_idx <- match(snp_IDs, names(tt))
  ids <- rep(snp_IDs, times = tt[tt_idx])

  idx <- match(ids, df$SNP)
  if (any(is.na(idx))) stop('bad SNP ids')
#  browser()
  res <- df[idx, ]
  rownames(res) <- NULL

  res
}


