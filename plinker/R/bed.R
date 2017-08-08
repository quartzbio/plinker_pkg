
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
#' @param bed			the .bed filename
#' @param bim			the .bim
#' @param fam 		the .fam path
#' @param fam_df	the fam data as a data frame
#' @param nb_snps	the number of snps
#' @keywords internal
new_bed <- function(bed, bim, fam, fam_df, nb_snps) {

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

  # the snp subset
  obj$snp_idx <- NULL

  class(obj) <- 'plinker_bed'


  obj
}


#' open a plink bed dataset
#'
#' @param prefix	the prefix path of the dataset (without the extensions)
#' @param files		the files of the dataset as a named list
#' @inheritParams new_bed
#' @return a plinker_bed object
#'
#' @export
bed_open <- function(prefix, files = as.list(unprefix_bed(prefix)),
  bed = files$bed, fam = files$fam, bim = files$bim)
{
  check_bed_files(bed, bim, fam)

  fam_df <- read_fam(fam)

  nb_snps <- infer_nb_snps(bed, nrow(fam_df))

  new_bed(bed, bim, fam, fam_df = fam_df, nb_snps = nb_snps)
}

### methods
#' @export
print.plinker_bed <- function(x, ...) {
  nb_snps1 <- bed_nb_snps(x, TRUE)
  nb_snps <- bed_nb_snps(x, FALSE)
  if (nb_snps != nb_snps1)
    nb_snps <- sprintf('%i/%i', nb_snps1, nb_snps)
  line <- sprintf('PLINK BED dataset: %s SNPs x %i samples\n',
    nb_snps, nrow(x$fam_df))
  cat(line)

  loaded <- if (exists('bim_df', envir = x$cache)) 'Loaded' else
          'Not yet loaded'
  line <- sprintf('bim annotations: "%s"\n', loaded)
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

#' get the number of samples in a plink bed dataset
#'
#' @inheritParams params
#' @family accessors
#' @export
bed_nb_samples <- function(bo) nrow(bo$fam_df)


#' get fam sample annotations
#'
#' @inheritParams params
#' @return the annotations as a data frame
#' @family accessors
#' @export
bed_fam_df <- function(bo) bo$fam_df

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
    assign('bim_df', bim_df, envir = bo$cache)
  }

  if (!subset || is.null(bo$snp_idx))
    bim_df
  else
    bim_df[bo$snp_idx, , drop = FALSE]
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
