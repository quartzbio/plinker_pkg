
#' read and parse a plink .fam file
#'
#' cf \url{http://www.cog-genomics.org/plink/1.9/formats#fam}
#'
#' @param path		the .fam file path
#' @return a data frame
#' @export
read_fam <- function(path) {
  utils::read.table(path,
    col.names = c('FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'SEX', 'PHENO'),
    colClasses = c(rep('character', 4), rep('integer', 2)))
}



#' save a plink .fam file
#'
#' cf \url{http://www.cog-genomics.org/plink/1.9/formats#fam}
#'
#' @param	df			the fam data frame to save
#' @param path		the path of the file to save
#' @param header	whether to write the header
#' @export
save_fam <- function(df, path, header = FALSE) {
  utils::write.table(df,path, sep = ' ', col.names = header, row.names = FALSE,
    quote = FALSE)
}


#' compute the sample IDs
#'
#' @param fam_df	the .fam data frame
#' @param ignore_fid	if TRUE the ids are the fam IID, otherwise it is
#' 	a concatenation of the fam FID and IID, separated by underscore
#' @return the ids
#' @family fam
#' @export
compute_sample_IDs <- function(fam_df, ignore_fid) {
  ids <- if (ignore_fid) {
      fam_df$IID
    } else {
      paste0(fam_df$FID, '_', fam_df$IID)
    }

  ids
}


#' merge a data frame with a FAM data frame
#'
#' a "left-join" style merge is performed on the FAM data frame.
#' The merge is performed:
#'   - either by (FID,IID) if present in df
#'   - or by sample_IDs versus id_var. in that case the id_var column
#'     is removed from the merged data frame
#'
#' N.B: the sample_IDs are computed from fam_df using ignore_fid
#'
#' will die unless each row of fam_df is matched exactly to one row of df
#'
#' columns from df with same names than columns from fam_df will be renamed with
#' a .y suffix
#'
#' the merged data frame is ordered as the original fam_df
#'
#' @param fam_df			the .fam data.frame, as returned by [read_fam()]. Not checked.
#' @param df					the data frame to merge
#' @param id_var			the df variable to use for the merge if (FID,IID) is not present
#' 										in df
#' @param ignore_fid	whether to ignore the .fam FID variable to compute the fam IDs
#' 										(cf [compute_sample_IDs()])
#' @return the merged data frame
#' @export
merge_df_with_fam <- function(fam_df, df, id_var = 'SUBJID', ignore_fid = FALSE) {
  cols <- c('FID', 'IID')
  to_remove <- NULL
  by.x <- by.y <- NULL

  if (all(cols %in% names(df))) {
    by.x <- by.y <- cols
  } else {
    if (!id_var %in% names(df)) stop('error, bad "id_var"')

    fam_df[[id_var]] <- compute_sample_IDs(fam_df, ignore_fid)

    to_remove <- by.x <- by.y <- id_var
  }

  res <- merge(fam_df, df, by.x = by.x, by.y = by.y,
    suffixes = c('', '.y'), sort = FALSE)

  if (!is.null(to_remove)) res[[to_remove]] <- NULL

  if (nrow(res) != nrow(fam_df))
    stop('bad merge, different number of rows in the output')

  res
}
