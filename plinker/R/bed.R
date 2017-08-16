
check_bed_files <- function(bed, bim, fam) {
  files <- list(bed = bed, bim = bim, fam = fam)
  for (i in seq_along(files)) {
    if (!file.exists(files[[i]])) {
      msg <- sprintf('Error, %s file "%s" does not exist.', names(files)[i],
        files[[i]])
      stop(msg)
    }
  }

  invisible()
}

#' create a bed object
#'
#' @param bed				the .bed filename
#' @param bim				the .bim
#' @param fam 			the .fam path
#' @param fam_df		the fam data as a data frame
#' @param nb_snps		the number of snps
#' @param ...				additional attributes to store in the object
#' @keywords internal
new_bed <- function(bed, bim, fam, fam_df, nb_snps, ...) {

  # user friendly check
  check_bed_files(bed, bim, fam)

  files <- list(bed = bed, bim = bim, fam = fam)
  obj <- normalizePath(unlist(files))
  names(obj) <- names(files)
  obj <- as.list(obj)

  obj$fam_df <- fam_df
  obj$nb_snps <- nb_snps


  # create the cache
  obj$cache <- new.env(parent = emptyenv())

  # add ...
  obj <- c(obj, list(...))

  class(obj) <- 'plinker_bed'


  obj
}


#' open a plink bed dataset
#'
#' @param prefix			the prefix path of the dataset (without the extensions)
#' @param files				the files of the dataset as a named list
#' @param ignore_fid	whether to ignore Family IDs in the .fam file. When used
#' 										in a non-family context, those ids are often dummy values
#' 										and are quite inconvenient. if NA, the default,
#' 										it will be TRUE iff the FIDs are constant/unique.
#' @inheritParams new_bed
#' @return a plinker_bed object
#'
#' @export
bed_open <- function(prefix, files = as.list(unprefix_bed(prefix)),
  bed = files$bed, fam = files$fam, bim = files$bim,
  ignore_fid = NA)
{
  check_bed_files(bed, bim, fam)

  fam_df <- read_fam(fam)

  if (is.na(ignore_fid)) {
    ignore_fid <- anyDuplicated(fam_df$FID) == 0
  }

  nb_snps <- infer_nb_snps(bed, nrow(fam_df))

  new_bed(bed, bim, fam, fam_df = fam_df, nb_snps = nb_snps,
    ignore_fid = ignore_fid)
}

### methods
#' @export
print.plinker_bed <- function(x, ...) {
  nb_snps1 <- bed_nb_snps(x, TRUE)
  nb_snps <- bed_nb_snps(x, FALSE)
  if (nb_snps != nb_snps1)
    nb_snps <- sprintf('%i/%i', nb_snps1, nb_snps)

  nbs1 <- bed_nb_samples(x, TRUE)
  nbs <- bed_nb_samples(x, FALSE)
  if (nbs != nbs1)
    nbs <- sprintf('%i/%i', nbs1, nbs)

  line <- sprintf('PLINK BED dataset: %s SNPs x %s samples\n', nb_snps, nbs)
  cat(line)

  loaded <- if (exists('bim_df', envir = x$cache)) 'Loaded' else
          'Not yet loaded'
  line <- sprintf('bim annotations: "%s"\n', loaded)
  cat(line)

  # ignore FIDs
  sids <- if (bed_ignore_fid(x)) 'IID' else 'FID_IID'
  line <- paste0('sample IDs: ', sids, '\n')
  cat(line)

  cat('---\n')
  paths <- x[c('bed', 'bim', 'fam')]
  cat(sprintf('%s:%s', names(paths), paths), sep = '\n')

  invisible()
}



### accessors
#' get the number of snps in a plink bed dataset
#'
#' @inheritParams params
#' @family accessors
#' @export
bed_nb_snps <- function(bo, subset = TRUE) {
  if (!subset || is.null(bo$snp_idx)) {
    bo$nb_snps
  } else {
    length(bo$snp_idx)
  }
}

#' get the current subset of snps as indices
#'
#' @inheritParams params
#' @return the indices, or NULL if not subset
#' @family accessors
#' @export
bed_snp_idx <- function(bo) bo$snp_idx

#' get the current subset of samples as indices
#'
#' @inheritParams params
#' @return the indices, or NULL if not subset
#' @family accessors
#' @export
bed_sample_idx <- function(bo) bo$sample_idx

#' get the current subset of samples as ids
#'
#' @inheritParams params
#' @param ignore_fid	if TRUE the ids are the fam IID, otherwise it is
#' 	a concatenation of the fam FID and IID, separated by underscore
#' @return the ids
#' @family accessors
#' @export
bed_sample_IDs <- function(bo, subset = TRUE, ignore_fid = bed_ignore_fid(bo)) {
  df <- bed_fam_df(bo, subset = subset)
  ids <- if (ignore_fid) {
    df$IID
  } else {
    paste0(df$FID, '_', df$IID)
  }

  ids
}


#' get the number of samples in a plink bed dataset
#'
#' @inheritParams params
#' @family accessors
#' @export
bed_nb_samples <- function(bo, subset = TRUE) {
  if (!subset || is.null(bo$sample_idx)) {
    nrow(bo$fam_df)
  } else {
    length(bo$sample_idx)
  }
}


#' whether the fam FIDs are ignored
#'
#' @inheritParams params
#' @return logical(1)
#' @family accessors
#' @export
bed_ignore_fid <- function(bo) { bo$ignore_fid }


#' get fam sample annotations
#'
#' @inheritParams params
#' @return the annotations as a data frame
#' @family accessors
#' @export
bed_fam_df <- function(bo, subset = TRUE) {
  if (!subset || is.null(bo$sample_idx))
    bo$fam_df
  else
    bo$fam_df[bo$sample_idx, , drop = FALSE]
}

#' get/fetch bim snp annotations lazyly
#'
#' For big datasets, reading the .bim annotations can take a long time.
#' The first time this function is called, it will actually read the file,
#' but all subsequent calls will use the cached value.
#'
#' @inheritParams params
#' @return the annotations as a data frame
#' @family accessors
#' @export
bed_bim_df <- function(bo, subset = TRUE) {
  bim_df <- get0('bim_df', envir = bo$cache, inherits = FALSE)
  if (is.null(bim_df)) {
    bim_df <- read_bim(bo$bim)
    if (nrow(bim_df) != bed_nb_snps(bo, subset = FALSE))
      stop('Error, the nb of SNPs in .bim is different from the inferred nb!')
    assign('bim_df', bim_df, envir = bo$cache)
  }

  if (!subset || is.null(bo$snp_idx))
    bim_df
  else
    bim_df[bo$snp_idx, , drop = FALSE]
}

#' subset  a plink dataset
#'
#' @inheritParams params
#' @return the subsetted bed dataset object
#' @family subset
#' @export
bed_subset <- function(bo,
  snp_IDs = NULL, snp_idx = NULL, sample_IDs = NULL, sample_idx = NULL)
{
  if (!is.null(snp_IDs) && !is.null(snp_idx))
    stop('incompatible parameters snp_IDs and snp_idx')

  if (!is.null(sample_IDs) && !is.null(sample_idx))
    stop('incompatible parameters sample_IDs and sample_idx')

  bo2 <- if (!is.null(snp_IDs)) {
    bed_subset_snps_by_IDs(bo, snp_IDs)
  } else if (!is.null(snp_idx)) {
    bed_subset_snps_by_idx(bo, snp_idx)
  } else { NULL }

  if (!is.null(bo2)) bo <- bo2

  bo2 <- if (!is.null(sample_IDs)) {
      bed_subset_samples_by_IDs(bo, sample_IDs)
  } else if (!is.null(sample_idx)) {
      bed_subset_samples_by_idx(bo, sample_idx)
  } else { NULL }

  if (!is.null(bo2)) bo <- bo2

  bo
}

infer_nb_snps <- function(bed, nb_samples) {
  (file.size(bed)  - 3L) / ceiling(nb_samples / 4L)
}



unprefix_bed <- function(prefix) {
  types <- c('bed', 'bim', 'fam')
  fns <- paste0(prefix, '.', types)
  names(fns) <- types

  fns
}



