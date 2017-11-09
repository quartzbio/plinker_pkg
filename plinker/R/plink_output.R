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



#' reorder if needed a sample-based plink output
#'
#' useful to apply the ordering of a bed_plink object
#'
#' @param df		a PLINK output as a data frame, that has FID,IID columns
#' @inheritParams params
#' @return the reordered output as a data frame
#' @keywords internal
reorder_plink_sample_output <- function(df, sample_IDs, ignore_fid = FALSE) {

  if (!all(c('FID', 'IID') %in% names(df))) stop('missing FID or IID column(s)')

  if (length(sample_IDs) == 0) return(df)

  ids <- compute_sample_IDs(df, ignore_fid)
  if (!all(ids %in% sample_IDs))  stop('bad sample ids')

  idx <- match(ids, sample_IDs)
  ord <- order(idx)

  res <- df[ord, ]
  rownames(res) <- NULL

  res
}

#' substitute the FAM ids (FID, IID) by other vars in a PLINK-like sample output
#'
#' @param df			a PLINK output as a data frame, that has FID,IID columns
#' @param annot	  a data frame with FID,IID and "id" columns
#' @inheritParams params
#' @return a merged data frame without the FID/IID columns, but the "id" column
#' 	instead
#' @keywords internal
annotate_plink_sample_output <- function(df, annot, id) {
  if (length(id) != 1 || ! id %in% names(annot))
    stop('bad param id')

  df2 <- merge_df_with_fam(df, annot)

  # reorder cols
  extra_cols <- setdiff(names(annot), id)
  newcols <- c(id, setdiff(names(df2), names(annot)))
  df2 <- df2[, newcols, drop = FALSE]

  df2
}

#' add the a annotaion id to a PLINK-like snp output
#'
#' @param df				a PLINK output as a data frame, that has a SNP column
#' @param annot	  	a data frame with a SNP and "id" columns
#' @param annot_id	the id var name in annot
#' @inheritParams params
#' @return a data frame with an additional "id" column next to the SNP col
#' @keywords internal
annotate_plink_snp_output <- function(df, annot, annot_id,
  annot_id_replacement = annot_id, df_snp_var = 'SNP') {
  if (length(annot_id) != 1 || ! annot_id %in% names(annot))
    stop('bad param annot_id')
  SNP <- 'SNP'
  cols <- c(SNP, annot_id)
  if (any(! cols %in% names(annot)))
    stop('mandatory cols not found in annot', paste0(cols, collapse = ', '))

  idx <- match(df[[df_snp_var]], annot$SNP)

  allcols <- names(df)
  df[[annot_id_replacement]] <- annot[[annot_id]][idx]

  # reorder cols
  newcols <- append(allcols, annot_id_replacement, match(df_snp_var, allcols))
  df <- df[, newcols, drop = FALSE]

  df
}


#' reorder if needed a SNP-based plink output
#'
#' useful to apply the ordering of a bed_plink object
#'
#' @param df		a PLINK output as a data frame, that has a SNP column
#' @inheritParams params
#' @return the reordered output as a data frame
#' @keywords internal
reorder_plink_snp_output <- function(df, snp_IDs) {
  if (!'SNP' %in% names(df)) stop('no SNP column')

  if (length(snp_IDs) == 0) return(df)

  if (!all(df$SNP %in% snp_IDs))  stop('bad SNP ids')

  idx <- match(df$SNP, snp_IDs)
  ord <- order(idx)

  res <- df[ord, ]
  rownames(res) <- NULL

  res
}


