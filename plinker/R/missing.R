
compute_missing_value_candidate_for_numeric <- function(v, candidate = -9L) {
  if (!is.numeric(v)) stop('v must be numeric')
  candidate <- as.integer(candidate)

  # always try -9L for legacy reason
  if (!any(v == candidate, na.rm = TRUE)) return(candidate)

  candidate <- as.integer(floor(min(v) - 0.5))
  if (candidate >= 0) candidate <- -9L

  candidate
}


compute_missing_value_candidate <- function(pheno, covars = NULL, ...) {
  values <- if (is.null(covars))
      pheno else
      c(pheno, as.numeric(data.matrix(covars)))
  compute_missing_value_candidate_for_numeric(values, ...)
}




