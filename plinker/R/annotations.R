
#' add custom sample annotations
#'
#' The annotations will be merged with the FAM data using [merge_df_with_fam()].
#' The merge_id is only used for the merge. The "id", if given, will be used
#' in all relevant output.
#'
#' @inheritParams params
#' @inheritParams merge_df_with_fam
#' @inheritDotParams merge_df_with_fam -fam_df -df
#' @param df		the sample annotations
#' @param id		an optional name of a column in df that is to be used in further
#' 		output as main custom sample ID. Its values must be character and unique
#' @return the annotated plinker_bed object
#' @export
bed_set_sample_annot <- function(bo, df, id = NULL, ...) {
  fam <- bed_fam_df(bo)
  annot <- merge_df_with_fam(fam, df, ignore_fid = bed_ignore_fid(bo), ...)

  cols_to_keep <- setdiff(names(annot), names(fam))
  annot <- annot[, cols_to_keep, drop = FALSE]

  bo$sample_annot <- annot

  if (!is.null(id)) {
    values <- df[[id]]
    if (is.null(values)) stop('bad column name id:', id)
    if (!is.character(values)) stop('"id" column must be a character vector:', id)
    if (anyDuplicated(values) > 0) stop('duplicated values in column :', id)

    bo$sample_annot_id <- id
  }


  bo
}

#' get custom sample annotations if any
#'
#' @inheritParams params
#' @param var		an optional variable name
#' @param fam		whether to include the FAM data (only if no var has been given
#' @return the annotation data frame, or just a column as a vector
#' 	if a var is given
#' @export
bed_get_sample_annot <- function(bo, var = NULL, fam = FALSE) {
  if (!is.null(var))
    return(bo$sample_annot[[var]])

  annot <- bo$sample_annot

  if (fam)
    annot <- cbind(bed_fam_df(bo), annot, stringsAsFactors = FALSE)

  annot
}

#' get the sample annotations (custom) ID if any
#'
#' @inheritParams params
#' @return the name of the ID var, or NULL if none
#' @export
bed_get_sample_annot_id <- function(bo) {
  bo$sample_annot_id
}


#' add custom snp annotations
#'
#' The annotations will be merged with the BIM data.
#'
#' @inheritParams params
#' @param df		the snp annotations, with a SNP column, and the same number of
#' 		rows than the BIM df.
#' @param id		an optional name of a column in df that is to be used in further
#' 		output as main custom SNP ID. Its values must be character.
#' @return the annotated plinker_bed object
#' @export
bed_set_snp_annot <- function(bo, df, id = NULL) {
  if (is.null(df$SNP)) stop('no SNP column')
  if (anyDuplicated(df$SNP) > 0) stop('duplicated values in SNP column')

  bim_df <- bed_bim_df(bo)

  res <- merge(bim_df[, 'SNP', drop = FALSE], df,
     by = 'SNP', suffixes = c('', '.y'), sort = FALSE)
  if (nrow(res) < nrow(bim_df)) stop('non exhaustive annotations')

  bo$snp_annot <- res

   if (!is.null(id)) {
      values <- res[[id]]
      if (is.null(values)) stop('bad column name id:', id)
      if (!is.character(values)) stop('"id" column must be a character vector:', id)
      if (anyDuplicated(values) > 0) warning('duplicated values in column :', id)

      bo$snp_annot_id <- id
  }

  bo
}

#' get custom snp annotations if any
#'
#' @inheritParams params
#' @return the annotation data frame
#' @export
bed_get_snp_annot <- function(bo) {
  bo$snp_annot
}

#' get the custom snp annotations ID if any
#'
#' @inheritParams params
#' @return the name of the ID var, or NULL if none
#' @export
bed_get_snp_annot_id <- function(bo) {
  bo$snp_annot_id
}

