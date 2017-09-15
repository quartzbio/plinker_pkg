#' check the covariates to be used with PLINK
#'
#' N.B: the missing value should (still) be encoded as NA
#'
#' @inheritParams params
#' @inheritParams merge_df_with_fam
#' @param covars 	the covars as a data frame,
#' 		as e.g. returned by [bed_make_covars()]. must not contain constant covars,
#' 		nor categorical covars.
#' @keywords internal
check_plink_covars <- function(covars, fam_df) {

  if (!identical(names(covars)[1:2], c('FID', 'IID'))) {
    stop('missing required first columns FID and IID')
  }
  # check IDs
  if (!all(mapply(identical, covars[, 1:2], fam_df[, 1:2]))) {
    stop('columns FID and IID do not match')
  }

  ### check covars columns content (TODO)

  # check for constant covars
  df <- covars[, -(1:2), drop = FALSE]
  .uniq_len <- function(col) length(stats::na.omit(unique(col)))
  lens <- sapply(df, .uniq_len)
  ones <- which(lens == 1L)
  if (length(ones) > 0) {
    stop('Error, constant covariates: ', paste0(names(df)[ones], collapse = ', '))
  }

  # check for categorical variables
  .is_categorical <- function(col) is.character(col) || is.factor(col)
  categs <- which(sapply(df, .is_categorical))
  if (length(categs) > 0) {
    stop('Error, categorical covariates: ',
      paste0(names(df)[categs], collapse = ', '))
  }

  invisible()
}

#' use model.matrix to dummy code categorical variables
#'
#' @param df	the data frame to recode
#'
#' @return a numeric matrix
#' @keywords internal
dummify_df_vars <- function(df, ...) {
  vars <- names(df)
  fo <- paste0('~', paste0(vars, collapse = ' + '))

  res <- stats::model.matrix(stats::as.formula(fo), df, ...)

  res[, -1, drop = FALSE]
}

#' format,recode and check covariates to fit the PLINK format
#'
#' PLINK expects covariates in a table with FID, IID and covariates columns
#' cf <https://www.cog-genomics.org/plink/1.9/input#covar>
#'
#' This function helps building this table
#'
#' @inheritParams params
#' @param covars 	the covars as a data frame, with also sufficient ID columns
#' 	to be merged with the .fam table (cf [merge_df_with_fam()])
#' @param ... passed to [stats::model.matrix()], so that you can customize
#' 	e.g. your contrasts
#'
#' @return the covariates formatted for PLINK, as a data frame
#' @export
bed_make_covars <- function(bo, covars, ...) {
  fam_df <- bed_fam_df(bo)
  df <- merge_df_with_fam(fam_df, covars, ignore_fid = bed_ignore_fid(bo), ...)

  df <- df[, -(3:6)]

  mat <- dummify_df_vars(df[, -(1:2), drop = FALSE], ...)

  df <- data.frame(df[, 1:2], as.data.frame(mat),
    stringsAsFactors = FALSE, row.names = NULL)

  df
}


