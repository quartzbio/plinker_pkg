#' check the covariates to be used with PLINK
#'
#' @inheritParams params
#' @inheritParams merge_df_with_fam
#' @param covars 	the covars as a data frame,
#' 		as e.g. returned by [bed_make_covars()]
#'
#' @return the covariates formatted for PLINK, as a data frame
#' @export
check_plink_covars <- function(covars, fam_df) {

  if (!identical(names(covars)[1:2], c('FID', 'IID'))) {
    stop('missing required first columns FID and IID')
  }
  # check IDs
  if (!identical(covars[, 1:2], fam_df[, 1:2])) {
    stop('columns FID and IID do not match')
  }

  # check covars columns content (TODO)


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

  res <- stats::model.matrix(as.formula(fo), df, ...)

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


