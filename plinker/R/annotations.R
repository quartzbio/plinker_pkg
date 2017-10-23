
#' add sample annotations
#'
#' The annotations will be merged with the FAM data using [merge_df_with_fam()].
#' The merge_id is only used for the merge. The real_id, if given, will be used
#' in all relevant output.
#'
#' @inheritParams params
#' @inheritParams merge_df_with_fam
#' @inheritDotParams merge_df_with_fam -fam_df -df
#' @param df		the sample annotations
#' @param id		an optional name of a column in df that is to be used in further
#' 		output as main secondary sample ID. Its values must be character and unique
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

#' get sample annotations if any
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

