


#' run a plink command on a plinker_bed object
#'
#' @param allow_no_sex		if FALSE, samples with ambiguous sex have their
#' 		 phenotypes set to missing, cf the \code{--allow-no-sex} option in
#' 			\url{http://www.cog-genomics.org/plink/1.9/filter#maf}
#' @param nonfounders		by default in plink, nonfounders are not counted
#' 		by --freq{x} or --maf/--max-maf/--hwe. setting this to TRUE disable
#' 		this behaviour by using the \code{--nonfounders} option, cf
#' 		\url{http://www.cog-genomics.org/plink/1.9/filter#nonfounders}
#'
#' @param keep_allele_order		whether to tell plink not to reorder alleles based
#' 								on maf for instance. if TRUE adds the --keep-allele-order
#' 								otherwise does nothing
#' @param lexicographic_allele_order
#' 														tell PLINK to use the lexicographic order for alleles
#' @param a2_alleles					tell plink to use the supplied alleles as A2 allele
#' @param missing_phenotype		tell PLINK which integer value is used to encode
#' 														missing phenotypes.
#' @param threads							max number of threads for PLINK to use.
#' 														N.B: not all algorithms are parallelized. if NA
#' 														do not use the --threads option. i.e. let PLINK
#' 														decide, cf <https://www.cog-genomics.org/plink/1.9/other#memory>
#'
#' @param ...		passed to \code{\link{plink_cmd}}
#' @inheritParams params
#' @inheritParams plink_cmd
#' @export
bed_plink_cmd <- function(bo,
  args,
  snp_idx = bed_snp_idx(bo),
  sample_idx = bed_sample_idx(bo),
  allow_no_sex = TRUE,
  nonfounders = TRUE,
  keep_allele_order = !lexicographic_allele_order,
  lexicographic_allele_order = FALSE,
  a2_alleles = if (lexicographic_allele_order) bed_allele_higher(bo) else NULL,
  missing_phenotype = -9L,
  threads = NA,
  ...)
{
  ### add plink input files
  fns <- bo[c('bed', 'fam', 'bim')]
  fargs <- paste0('--', names(fns), ' ', fns)
  args <- c(fargs, args)

  if (!is.null(a2_alleles)) {
    if (length(a2_alleles) != bed_nb_snps(bo))
      stop('bad a2_alleles length', length(a2_alleles))
    a2_fn <- 'input_a2_alleles.txt'
    bim2 <- bed_bim_df(bo)
    bim2$NEWA2 <- a2_alleles
    save_bim(bim2[, c('SNP', 'NEWA2')], a2_fn)
    args <- c(paste('--a2-allele', a2_fn), args)
  }

  ## add optional flags
  if (allow_no_sex) args <- c(args, '--allow-no-sex')
  if (nonfounders) args <- c(args, '--nonfounders')
  if (keep_allele_order) args <- c(args, '--keep-allele-order')
  if (!is.na(threads)) args <- c(args, paste0('--threads ', threads))

  check_missing_phenotype(missing_phenotype)
  args <- c(args, sprintf('--missing-phenotype %i', missing_phenotype))

  ### add subsetting args if needed
  plink_snp_ids_fn <- 'input_snp_ids.txt'
  if (!is.null(snp_idx)) {
    bim_df <- bed_bim_df(bo, subset = FALSE)
    snp_ids <- bim_df$SNP[snp_idx]

    writeLines(snp_ids, plink_snp_ids_fn)

    args <- c(paste('--extract', plink_snp_ids_fn), args)
  }

  plink_sample_ids_fn <- 'input_sample_ids.txt'
  if (!is.null(sample_idx)) {
    fam_df <- bed_fam_df(bo, subset = FALSE)

    readr::write_tsv(fam_df[sample_idx, 1:2, drop = FALSE],
      plink_sample_ids_fn, col_names = FALSE)

    args <- c(paste('--keep', plink_sample_ids_fn), args)
  }

  invisible(plink_cmd(args, ...))
}

#' compute genotypic distances between samples
#'
#' N.B: the default PLINK missingness correction is disabled for now
#' because I was not able to reproduce it.
#'
#' Note: is is multithreaded, but not very efficiently. In my benchmarks, after
#' two threads performance degrades very badly.
#'
#' @param type 	distance unit type.
#' 	cf <https://www.cog-genomics.org/plink/1.9/distance#read_dists>
#' @param missingness		the missingness correction to apply for missing genotypes
#' 	Currently only flat-missing
#' @param ...		passed to \code{\link{bed_plink_cmd}}
#' @inheritParams bed_plink_cmd
#' @return a distance matrix
#'
#' @seealso bed_plink_cmd
#' @export
bed_plink_distance <- function(bo, type = c('ibs', 'allele-ct'),
  missingness = c('flat-missing'), threads = 1, ...)
{
  setup_temp_dir()

  type <- match.arg(type)
  missingness <- match.arg(missingness)

  args <- c('--distance square bin', type)
  args <- c(args, missingness)
  bed_plink_cmd(bo, args, threads = threads, ...)

  output <- switch(type, ibs = 'mibs', 'dist')
  id_file <- sprintf('plink.%s.id', output)
  dist_file <- sprintf('plink.%s.bin', output)

  df <- read_plink_output(id_file, header = FALSE)
  nb <- nrow(df)

  dist <- readBin(dist_file, what = 'numeric', n = nb * nb)
  dim(dist) <- c(nb, nb)

  ### reorder matrix
  ids <- paste0(df[[1]], '_', df[[2]])
  ordered_ids <- bed_sample_IDs(bo, ignore_fid = FALSE)

  idx <- match(ordered_ids, ids)
  dist2 <- dist[idx, idx]

  rownames(dist2) <- bed_sample_IDs(bo, custom = TRUE)
  colnames(dist2) <- rownames(dist2)

  dist2
}


#' compute allele frequencies using plink --freq counts
#'
#' N.B: plink --freq does NOT reorder the alleles.
#'
#' @param ...		passed to \code{\link{bed_plink_cmd}}
#' @inheritParams bed_plink_cmd
#' @return a data frame, cf \url{http://www.cog-genomics.org/plink/1.9/formats#frq_count}
#'
#' @seealso bed_plink_cmd
#' @export
bed_plink_freq_count <- function(bo, ...)
 {
  setup_temp_dir()

  args <- '--freq counts'
  bed_plink_cmd(bo, args, ...)
  df <- read_plink_freq_counts('plink.frq.counts')

  postprocess_snp_output(bo, df)
}

#' compute LD using plink --r2 dprime
#'
#' N.B: use min_r2 to limit output !!!!
#'
#' Currently the implementation does not seem parallelized.
#'
#' This does not depend on keep_allele_order.
#'
#' @inheritDotParams bed_plink_cmd -bo -args
#' @inheritParams params
#' @param window_size		max number of SNPs in LD window
#' @param window_length	max window length in kb
#' @param min_r2				minimum r2 value to report. Very important for performance
#'                      and to avoid filling up your hard-disk.
#' @return a data frame with columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2 DP,
#' 	cf <https://www.cog-genomics.org/plink/1.9/ld>
#' @export
#' @md
bed_plink_ld <- function(
  bo,
  window_size = 10L,
  window_length = 1000000L,
  min_r2 = 0,
  ...
) {

  main_args <- '--r2 dprime'
  window_args <- sprintf('--ld-window %i --ld-window-kb %i --ld-window-r2 %s',
    window_size, window_length, min_r2)

  args <- paste(main_args, window_args, sep = ' ')

  dir <- setup_temp_dir()
  bed_plink_cmd(bo, args, ...)
  ld <- read_plink_output('plink.ld')

  if (!is.null(id <- bed_get_snp_annot_id(bo))) {
    ld <- annotate_plink_snp_output(ld, bed_get_snp_annot(bo), id,
      df_snp_var = 'SNP_A', annot_id_replacement = paste0(id, '_A'))
    ld <- annotate_plink_snp_output(ld, bed_get_snp_annot(bo), id,
      df_snp_var = 'SNP_B', annot_id_replacement = paste0(id, '_B'))
  }

  ld
}

#' select a subset of SNPS that are in approximate linkage equilibrium with each
#'  other using plink --indep-pairphase
#'
#' cf <https://www.cog-genomics.org/plink/1.9/ld#indep>
#'
#' Currently the implementation does not seem parallelized.
#'
#' This does not depend on keep_allele_order.
#'
#' @param step_size       the number of snps to slide the window by
#' @inheritParams bed_plink_ld
#' @inheritDotParams bed_plink_cmd -bo -args
#' @inheritParams params
#' @return the IDs of the selected sites as a character vector. N.B: these ids
#' 		have no particular order
#' @export
bed_plink_ld_select <- function(
  bo,
  window_size,
  window_length,
  step_size,
  min_r2,
  ...
) {

  # quick check args
  if (missing(window_size) && missing(window_length))
    stop('missing arg: window_size or window_length')

  if (!(missing(window_size) || missing(window_length)))
    stop('incompatible args: window_size and window_length')

  window_arg <- if (!missing(window_size)) as.character(window_size) else
      sprintf('%skb', window_length)

  # brand new option: --indep-pairphase
  args <- sprintf('--indep-pairphase %s %i %s', window_arg, step_size,
    min_r2)

  dir <- setup_temp_dir()
  bed_plink_cmd(bo, args, ...)

  ids <- readLines('plink.prune.in')

  ids
}

#' compute Fisher tests using PLINK --model fisher
#'
#' Tests performed:
#'   - 1df chi-square allelic test
#'   - 1df dominant gene action
#'   - 1df recessive gene action
#'   - 2df genotypic
#'   - Cochran-Armitage trend
#'
#' ** The results depend on keep_allele_order **.
#'
#' cf <https://www.cog-genomics.org/plink/1.9/assoc#model>
#'
#' Currently the implementation is not parallelized.
#'
#' @inheritDotParams bed_plink_cmd -bo -args
#' @inheritParams params
#' @param phenotype		a binary phenotype vector as an integer vector
#' 	Case/control phenotypes are expected to be encoded as:
#'   - 1=unaffected (control)
#'   - 2=affected (case)
#' @export
#' @family plink
#' @seealso bed_phenotype_from_df
bed_plink_fisher <- function(
  bo,
  phenotype = bed_fam_df(bo)$PHENO,
  ...
) {
  check_case_control_phenotype(phenotype, bed_nb_samples(bo))
  missing_phenotype <- compute_missing_value_candidate(phenotype)
  phenotype[is.na(phenotype)] <- missing_phenotype

  defined <- sum(phenotype %in% 1:2)
  if (defined == 0) stop('no non-missing phenotypes!')

  setup_temp_dir()

  # phenotype management
  pheno_path <- 'pheno.input'
  save_plink_phenotype_input(phenotype, bed_fam_df(bo), pheno_path)
  pheno_args <- sprintf('--pheno %s', pheno_path)

  args <- c('--model fisher', pheno_args)

  bed_plink_cmd(bo, args, missing_phenotype = missing_phenotype, ...)

  df <- read_plink_output('plink.model')

  postprocess_snp_output(bo, df)
}


#' compute linear regression using PLINK with covariates
#'
#'
#' Issues:
#'   - if a covar is constant PLINK can not compute anything. This function
#'    automatically detects this case and dies.
#'
#' Caution:
#'   - by default PLINK will reorder the alleles, and it changes the results.
#'     If you want to compare the results with another scan or subset you can use
#' 			the keep_allele_order param
#'
#' @inheritDotParams bed_plink_cmd -bo -args
#' @inheritParams params
#' @inheritParams recode_genotypes
#' @param phenotype a quantitative (unless logistic == TRUE)
#'        phenotype/response vector
#' @param covars	covariates as a data frame with first columns FID/IID
#' 	in the same order as bed_fam_df(bo). Use [bed_make_covars()] to easily
#' 	format them
#' @param logistic	whether to perform a logistic regression.
#' @export
#' @family plink
#' @seealso bed_phenotype_from_df
bed_plink_lm <- function(
  bo,
  phenotype = bed_fam_df(bo)$PHENO,
  model = c('', 'dominant', 'recessive'),
  covars = NULL,
  logistic = FALSE,
  quiet = FALSE,
  ...
) {
  info <- if (quiet) function(...) {}  else message
  setup_temp_dir()
  args <- sprintf('--linear %s', model)

  model <- match.arg(model)

  type <- if (logistic) 'logistic' else 'linear'
  args <- sprintf('--%s %s', type, model)

  ### covariates management:
  if (!is.null(covars)) {
    path <- 'covars.input'
    missing_phenotype <- compute_missing_value_candidate(phenotype,
      covars[, -(1:2), drop = FALSE])

    # apply it
    covars[is.na(covars)] <- missing_phenotype

    check_plink_covars(covars, bed_fam_df(bo))
    save_fam(covars, path, header = TRUE)
    covargs <- sprintf('--covar %s', path)
    args <- c(args, covargs)
  } else {
    missing_phenotype <- compute_missing_value_candidate(phenotype)
  }
  info('computed safe missing phenotype code: ', missing_phenotype)

  ### phenotype management:
  # N.B, mut be after covars because we need the final missing_phenotype value
  phenotype[is.na(phenotype)] <- missing_phenotype
  if (logistic)
    check_case_control_phenotype(phenotype, bed_nb_samples(bo))
  else
    check_quantitative_phenotype(phenotype, bed_nb_samples(bo))

  pheno_path <- 'pheno.input'
  save_plink_phenotype_input(phenotype, bed_fam_df(bo), pheno_path)
  pheno_args <- sprintf('--pheno %s', pheno_path)
  args <- c(args, pheno_args)


  bed_plink_cmd(bo, args, missing_phenotype = missing_phenotype, quiet = quiet,
    ...)

  df <- read_plink_output(sprintf('plink.assoc.%s', type))

  postprocess_snp_output(bo, df)
}


save_plink_phenotype_input <- function(phenotype, fam_df, path) {
  fam_df$PHENO <- phenotype
  ### with the --pheno option, PLINK reads the phenotype from the 3rd column
  pheno_df <- fam_df[, c('FID', 'IID', 'PHENO')]
  save_fam(pheno_df, path)
}


#' compute genotype frequencies using plink --freqx
#'
#' @inheritDotParams bed_plink_cmd -bo -args
#' @inheritParams params
#' @param keep_allele_order cf the ... param
#'
#' @return a data frame, cf <http://www.cog-genomics.org/plink/1.9/formats#frqx>
#'
#' @seealso bed_plink_cmd
#' @export
#' @md
bed_plink_freqx <- function(bo, keep_allele_order = TRUE, ...)
{
  setup_temp_dir()

  args <- '--freqx'
  bed_plink_cmd(bo, args, keep_allele_order = keep_allele_order, ...)
  df <- read_plink_output('plink.frqx')

  postprocess_snp_output(bo, df)
}

#' compute snp and sample missing rates using plink --missing
#'
#' @param ...		passed to \code{\link{bed_plink_cmd}}
#' @inheritParams bed_plink_cmd
#' @return a list of 2 data frame,
#' 	cf [http://www.cog-genomics.org/plink/1.9/formats#imiss](imiss) and
#' 	[http://www.cog-genomics.org/plink/1.9/formats#lmiss](lmiss)
#'
#' @seealso bed_plink_cmd
#' @export
#' @md
bed_plink_missing <- function(bo, ...)
{
  setup_temp_dir()

  args <- '--missing'
  bed_plink_cmd(bo, args, ...)
  imiss <- read_plink_output('plink.imiss')
  lmiss <- read_plink_output('plink.lmiss')

  lmiss <- postprocess_snp_output(bo, lmiss)

  imiss <- postprocess_sample_output(bo, imiss)

  list(sample = imiss, snp = lmiss)
}

postprocess_sample_output <- function(bo, out) {
  out <- reorder_plink_sample_output(out, bed_sample_IDs(bo), ignore_fid = bed_ignore_fid(bo))
  if (!is.null(id <- bed_get_sample_annot_id(bo))) {
    out <- annotate_plink_sample_output(out,
      bed_get_sample_annot(bo, fam = TRUE), id)
  }

  out
}

postprocess_snp_output <- function(bo, out, ...) {
  out <- reorder_plink_snp_output(out, bed_snp_IDs(bo))
  if (!is.null(id <- bed_get_snp_annot_id(bo))) {
    out <- annotate_plink_snp_output(out, bed_get_snp_annot(bo), id, ...)
  }

  out
}

#' output a ped file
#'
#' N.B: will also generate other files (.log, .map, .nosex)
#'
#' CAUTION: will not preserve the ordering of the subsets
#'
#' @inheritParams bed_plink_cmd
#' @param output_prefix	path of the .ped file to generate  (without the .ped suffix)
#' @param ...		passed to [bed_plink_cmd()]
#' @seealso bed_plink_cmd
#' @export
#' @md
bed_plink_ped <- function(bo, output_prefix, keep_allele_order = TRUE, ...)
{
  args <- c('--recode', paste0('--out ', output_prefix))
  bed_plink_cmd(bo, args, keep_allele_order = keep_allele_order, ...)
}


