context('bed_plink')


.bed_plink_distance <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  bo2 <- bed_subset(bo, sample_idx = 1:5)

  dist <- bed_plink_distance(bo2, quiet = TRUE)
  dist2 <- bed_R_distance(bo2)

  expect_equal(dist2, dist)

  dist <- bed_plink_distance(bo, quiet = TRUE)
  dist2 <- bed_R_distance(bo)
  expect_equal(dist2, dist)


  ### reordered subsets and sample annot
  ids <- bed_sample_IDs(bo)
  annot <- data.frame(
    MERGE_ID = rev(ids),
    SUBJID = paste0('ID_', rev(ids)),
    stringsAsFactors = FALSE)
  bo2 <- bed_set_sample_annot(bo, annot, id = 'SUBJID')
  bo2 <- bed_subset(bo2, sample_idx = 10:5, snp_idx = 10:6)

  dist <- bed_plink_distance(bo2, quiet = TRUE)
  dist2 <- bed_R_distance(bo2)
  expect_equal(dist2, dist)
  expect_identical(rownames(dist), bed_sample_IDs(bo2, custom = TRUE))
  expect_identical(colnames(dist), rownames(dist))


#  delta <- abs(dist2 - dist)
#  idx <- which.max(delta)
#
#  coords <- arrayInd(idx, dim(dist))
#  row <- coords[1, 1]
#  col <- coords[1, 2]
#  dist[row, col]
#  # [1] 0.5489331
#  dist2[row, col]
#  # [1] 0.53125
#
#  genos <- t(bed_genotypes(bo))
#
#  g1 <- unname(genos[, row])
#  g2 <- unname(genos[, col])
#
#  state <- plinker:::ibs_state(g1, g2)
#  #  [1]  2  1  1  1  1  0  1  1  1  1  1  1  2 NA  1  0  2
#  plinker:::ibs(state)
#  # [1] 0.53125
#  score <- sum(state, na.rm = TRUE)
#  # [1] 17
#  score / (2 * length(state))
#  # [1] 0.5
#  score / (2*17)
#  # [1] 0.5
#  score / (2*16) # mine
#  # [1] 0.53125
#
#  score / (2*15)
#  # [1] 0.5666667
#  score / 31
#  # [1] 0.5483871
#
#  score / 32
#  # [1] 0.53125
#
#
#  18/ (32 + 2)
#  # [1] 0.5294118
#
#  dist[row, col]*30
#  # [1] 16.46799
#  dist[row, col]*32
#  # [1] 17.56586
#
#  bo3 <- bed_subset(bo, sample_idx = c(row, col))
#  d <- bed_plink_distance(bo3, quiet = F)
#  d2 <- bed_R_distance(bo3)
#
#
#  dist[idx]
#  dist2[idx]
#
#  expect_equal(dist2, dist)
#

}
test_that('bed_plink_distance', .bed_plink_distance())



.bed_plink_ld_select <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)

  # check incompatible params
  expect_error(bed_plink_ld_select(bo, window_size = 10, window_length = 100,
      min_r2 = 0.8), 'incompatible args')
  expect_error(bed_plink_ld_select(bo, min_r2 = 0.8), 'missing arg')

  # extreme setting => select one per chromosome !
  ws <- 20000
  ids <- bed_plink_ld_select(bo, window_size = ws, step_size = round(ws / 2),
    min_r2 = 0, quiet = TRUE)
  expect_length(ids, 1)

  # opposite: select all
  ids <- bed_plink_ld_select(bo, window_size = ws, step_size = round(ws / 2),
    min_r2 = 0.99, quiet = TRUE)
  expect_length(ids, nsnp)

  res_step1 <- bed_plink_ld_select(bo, window_size = 2, step_size = 1,
    min_r2 = 0.5, quiet = TRUE)
  res2_step2 <- bed_plink_ld_select(bo, window_size = 2, step_size = 2,
    min_r2 = 0.5, quiet = TRUE)

  ## because of step_size == 2, some snps are missed
  expect_gt(length(res2_step2), length(res_step1))

  res1 <- bed_plink_ld_select(bo, window_size = 2, step_size = 1,
    min_r2 = 0.5, quiet = TRUE)
  res2 <- bed_plink_ld_select(bo, window_size = 2, step_size = 1,
    min_r2 = 0.5, quiet = TRUE, lexicographic_allele_order = TRUE)
  expect_identical(res2, res1)
}
test_that('bed_plink_ld_select', .bed_plink_ld_select())



.bed_plink_ld <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  nsnp <- bed_nb_snps(bo)

  ld1 <- bed_plink_ld(bo, window_size = 2, quiet = TRUE)
  expect_equal(nrow(ld1), nsnp - 1L)
  expect_true(all(ld1$R2 >= 0 & ld1$R2 <=1))
  expect_true(all(ld1$DP >= 0 & ld1$DP <=1))

  ld2 <- bed_plink_ld(bo, window_size = 2,  keep_allele_order = TRUE,
    quiet = TRUE)
  expect_identical(ld2, ld1)

  ld3 <- bed_plink_ld(bo, window_size = 2,  lexicographic_allele_order = TRUE,
    quiet = TRUE)
  expect_identical(ld3, ld1)

  ld <- bed_plink_ld(bo, window_size = 3, quiet = TRUE)
  expect_equal(nrow(ld), nsnp * 2 - 3)

  res <- merge(ld1, ld, by = c('SNP_A', 'SNP_B'))
  expect_equal(res$R2.x, res$R2.y)
  expect_equal(res$DP.x, res$DP.y)

  ### window size
  pos <- bed_bim_df(bo)$POS

  expect_gt(max(ld1$BP_B - ld1$BP_A), 40000)
  ld <- bed_plink_ld(bo, window_size = 1000, window_length = 5, quiet = TRUE)
  expect_lt(max(ld$BP_B - ld$BP_A), 5001)

  ### min_r2
  ld <- bed_plink_ld(bo, window_size = 2, min_r2 = 0.5, quiet = TRUE)
  expect_equal(nrow(ld), sum(ld1$R2 >= 0.5))

  ### snp annot
  bim <- bed_bim_df(bo)
  annot <- data.frame(
    SNP = bim$SNP,
    ID = paste0('ID_', bim$SNP),
    TOTO = 1,
    stringsAsFactors = FALSE)

  bo2 <- bed_set_snp_annot(bo, annot, 'ID')
  res <- bed_plink_ld(bo, window_size = 2, quiet = TRUE)
  res2 <- bed_plink_ld(bo2, window_size = 2, quiet = TRUE)

  expect_identical(res2[, -c(4,8)], res)
  expect_identical(paste0('ID_', res2$SNP_A), res2$ID_A)
  expect_identical(paste0('ID_', res2$SNP_B), res2$ID_B)
}
test_that('bed_plink_ld', .bed_plink_ld())



.bed_plink_ped <- function() {
  setup_temp_dir()
  bo <- bed_open(plinker:::fetch_sample_bed())

  path <- 'toto'
  bed_plink_ped(bo, path, quiet = TRUE)
  ped <- paste0(path, '.ped')
  expect_true(file.exists(ped))

  df <- read_plink_ped(ped)
  expect_equal(nrow(df), bed_nb_samples(bo))
  expect_equal(ncol(df), bed_nb_snps(bo)*2 + 6)

  ##########

  bo2 <- bed_subset(bo, snp_idx = 10:5, sample_idx = 40:31)

  bed_plink_ped(bo2, path, keep_allele_order = TRUE, quiet = TRUE)
  df1 <- read_plink_ped(ped)

  bed_plink_ped(bo2, path, keep_allele_order = FALSE, quiet = TRUE)
  df2 <- read_plink_ped(ped)

  strs1 <- plinker:::extract_genotypes_from_ped(df1)
  strs2 <- plinker:::extract_genotypes_from_ped(df2)

  expect_identical( strs1[3, 4], "G/A")
  expect_identical( strs2[3, 4], "A/G")


  mat <- bed_genotypes_as_strings(bo2, sort = FALSE)
  ### BEWARE: bed_plink_ped loses the ordering !!
  mat2 <- mat[nrow(mat):1, ncol(mat):1]
  expect_equivalent(mat2, strs1)

  ### allele order: permute A1 and A2 in the ped output
  bed_plink_ped(bo, path, quiet = TRUE)
  df1 <- read_plink_ped(ped)
  bed_plink_ped(bo, path, lexicographic_allele_order = TRUE, quiet = TRUE)
  df2 <- read_plink_ped(ped)

  expect_equivalent(df1[6, 9:10], c('G', 'C'))
  expect_equivalent(df2[6, 9:10], c('C', 'G'))
}
test_that('bed_plink_ped', .bed_plink_ped())



.bed_plink_missing <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())
  fam <- bed_fam_df(bo)
  bim <- bed_bim_df(bo)

  res <- bed_plink_missing(bo, quiet = TRUE)

  imiss <- res$sample
  expect_identical(imiss[, 1:2], fam[, 1:2])
  expect_equal(unique(imiss$N_GENO), 17)
  expect_equal(max(imiss$N_MISS), 1)
  expect_equal(max(imiss$F_MISS), 1/17, tolerance = 1e-5)

  lmiss <- res$snp

  expect_identical(lmiss[, 1:2], bim[, 1:2])
  gfreqs <- bed_plink_freqx(bo, quiet = TRUE)
  expect_equal(gfreqs[[10]], lmiss$N_MISS)
  expect_equal(gfreqs[[10]] / 89, lmiss$F_MISS, tolerance = 1e-5)

  ######### ordering
  # SNP ordering
  bo2 <- bed_subset(bo, snp_IDs = rev(bed_snp_IDs(bo)))
  res <- bed_plink_missing(bo2, quiet = TRUE)
  ref <- bed_plink_missing(bo, quiet = TRUE)

  # no sample reordering here
  expect_identical(res$sample, ref$sample)

  expect_false(identical(res$snp, ref$snp))
  snp2 <- plinker:::reorder_plink_snp_output(res$snp, bed_snp_IDs(bo))
  expect_identical(snp2, ref$snp)

  # sample reordering
  bo2 <- bed_subset(bo, sample_IDs = rev(bed_sample_IDs(bo)))
  res <- bed_plink_missing(bo2, quiet = TRUE)
  ref <- bed_plink_missing(bo, quiet = TRUE)

  expect_false(identical(res$sample, ref$sample))
  samp2 <- plinker:::reorder_plink_sample_output(res$sample, bed_sample_IDs(bo),
    bed_ignore_fid(bo))
  expect_identical(samp2, ref$sample)

  ### allele order
  res1 <- bed_plink_missing(bo, quiet = TRUE)
  res2 <- bed_plink_missing(bo, quiet = TRUE, lexicographic_allele_order = TRUE)

  expect_identical(res2, res1)

  ### custom sample annotations
  ids <- bed_sample_IDs(bo)
  annot <- data.frame(
    MERGE_ID = rev(ids),
    SUBJID = paste0('ID_', rev(ids)),
    stringsAsFactors = FALSE)
  bo2 <- bed_set_sample_annot(bo, annot, id = 'SUBJID')
  res <- bed_plink_missing(bo2, quiet = TRUE)$sample
  ref <- bed_plink_missing(bo, quiet = TRUE)$sample

  # the output without the ids
  expect_identical(res[, -1], ref[, -(1:2)])

  # check the ids
  expect_identical(res$SUBJID, paste0('ID_', ref$FID, '_', ref$IID))

  ### sample reordering + sample annotation
  bo2 <- bed_subset(bo, sample_IDs = rev(bed_sample_IDs(bo)))
  bo3 <- bed_set_sample_annot(bo2, annot, id = 'SUBJID')
  res <- bed_plink_missing(bo3, quiet = TRUE)$sample
  ref <- bed_plink_missing(bo2, quiet = TRUE)$sample

  expect_identical(res[, -1], ref[, -(1:2)])
  expect_identical(res$SUBJID, paste0('ID_', ref$FID, '_', ref$IID))

  ### snp annotations
  annot <- data.frame(
    SNP = bim$SNP,
    ID = paste0('ID_', bim$SNP),
    TOTO = 1,
    stringsAsFactors = FALSE)

  bo2 <- bed_set_snp_annot(bo, annot, 'ID')

  res <- bed_plink_missing(bo, quiet = TRUE)$snp
  res2 <- bed_plink_missing(bo2, quiet = TRUE)$snp

  expect_identical(res2[, -3], res)
  expect_identical(paste0('ID_', res2$SNP), res2$ID)
}
test_that('bed_plink_missing', .bed_plink_missing())




.bed_plink_freqx <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  bim_df <- bed_bim_df(bo)

  ### no subset
  df <- bed_plink_freqx(bo, quiet = TRUE)

  nb <- unique(rowSums(df[, 5:10]))
  expect_equal(nb, bed_nb_samples(bo))

  expect_identical(df$SNP, bim_df$SNP)
  expect_identical(df$A1, bim_df$A1)
  expect_identical(df$A2, bim_df$A2)

  # sample filtering
  df2 <- bed_plink_freqx(bo, allow_no_sex = FALSE, quiet = TRUE)
  # only ignore phenotype for no sex, does not impact the freqs
  expect_identical(df2, df)

  df2 <- bed_plink_freqx(bo, nonfounders = FALSE, quiet = TRUE)
  nb <- unique(rowSums(df2[, 5:10]))
  nb_nonfounders <- bed_nb_samples(bo) - nb
  expect_equal(nb_nonfounders, 3)

  ## !! 3 non founders, although only 2 actual ones since on hasa parent not
  ## in the fam file

  ### ordering
  bo2 <- bed_subset(bo, snp_IDs = rev(bed_snp_IDs(bo)))
  res <- bed_plink_freqx(bo2, quiet = TRUE)
  ref <- bed_plink_freqx(bo, quiet = TRUE)

  expect_false(identical(res, ref))
  res2 <- plinker:::reorder_plink_snp_output(res, bed_snp_IDs(bo))
  expect_identical(res2, ref)

  ### allele order
  df1 <- bed_plink_freqx(bo, quiet = TRUE)
  df2 <- bed_plink_freqx(bo, quiet = TRUE, lexicographic_allele_order = TRUE)

  expect_true(df1$A2[2] != df2$A2[2])

  expect_true(df1[2, 5] == df2[2, 7])
  expect_true(df1[2, 5] == df2[2, 7])
  expect_true(df1[2, 6] == df2[2, 6])

  ### snp annot
  bim <- bed_bim_df(bo)
  annot <- data.frame(
    SNP = bim$SNP,
    ID = paste0('ID_', bim$SNP),
    TOTO = 1,
    stringsAsFactors = FALSE)

  bo2 <- bed_set_snp_annot(bo, annot, 'ID')
  res <- bed_plink_freqx(bo, quiet = TRUE)
  res2 <- bed_plink_freqx(bo2, quiet = TRUE)

  expect_identical(res2[, -3], res)
  expect_identical(paste0('ID_', res2$SNP), res2$ID)

}
test_that('bed_plink_freqx', .bed_plink_freqx())



.bed_plink_freq_count <- function() {
  bo <- bed_open(plinker:::fetch_sample_bed())

  bim_df <- bed_bim_df(bo)

  ### no subset
  df <- bed_plink_freq_count(bo, quiet = TRUE)

  nb <- unique(with(df, C1 + C2 + 2*G0)) / 2
  expect_equal(nb, bed_nb_samples(bo))

  expect_identical(df$SNP, bim_df$SNP)
  expect_identical(df$A1, bim_df$A1)
  expect_identical(df$A2, bim_df$A2)

  # sample filtering
  df2 <- bed_plink_freq_count(bo, allow_no_sex = FALSE, quiet = TRUE)
  # only ignore phenotype for no sex, does not impact the freqs
  expect_identical(df2, df)

  df2 <- bed_plink_freq_count(bo, nonfounders = FALSE, quiet = TRUE)
  nb <- unique(with(df2, C1 + C2 + 2*G0)) / 2
  nb_nonfounders <- bed_nb_samples(bo) - nb
  expect_equal(nb_nonfounders, 3)

  ## !! 3 non founders, although only 2 actual ones since on hasa parent not
  ## in the fam file

  #### subset by snps
  bo2 <- bed_subset_snps_by_idx(bo, 6:10)
  df2 <- bed_plink_freq_count(bo2, quiet = TRUE)
  expect_equivalent(df2, df[6:10, ])

  ### subset by samples
  bo2 <- bed_subset_samples_by_idx(bo, 21:60)
  df2 <- bed_plink_freq_count(bo2, quiet = TRUE)
  nb <- unique(with(df2, C1 + C2 + 2*G0)) / 2
  expect_equal(nb, 40)

  ### ordering
  bo2 <- bed_subset(bo, snp_idx = 10:2)
  df2 <- bed_plink_freq_count(bo2, quiet = TRUE)
  expect_identical(df2$SNP, bed_snp_IDs(bo2))

  ### allele order
  df1 <- bed_plink_freq_count(bo, quiet = TRUE)
  df2 <- bed_plink_freq_count(bo, quiet = TRUE, lexicographic_allele_order = TRUE)

  expect_true(df1$A2[2] != df2$A2[2])

  ### BUG!!!!!!: PLINK --freq counts is buggy
  # regression test for bug now fixed in PLINK
  # expect_true(df1$C2[2] == df2$C2[2])

  ### snp annot
  bim <- bed_bim_df(bo)
  annot <- data.frame(
    SNP = bim$SNP,
    ID = paste0('ID_', bim$SNP),
    TOTO = 1,
    stringsAsFactors = FALSE)

  bo2 <- bed_set_snp_annot(bo, annot, 'ID')
  res <- bed_plink_freq_count(bo, quiet = TRUE)
  res2 <- bed_plink_freq_count(bo2, quiet = TRUE)

  expect_identical(res2[, -3], res)
  expect_identical(paste0('ID_', res2$SNP), res2$ID)
}
test_that('bed_plink_freq_count', .bed_plink_freq_count())



.bed_plink_cmd <- function() {
  plink_cmd <- plinker:::plink_cmd
  read_plink_freq <- plinker:::read_plink_freq

  bo <- bed_open(plinker:::fetch_sample_bed())

  setup_temp_dir()
  out <- bed_plink_cmd(bo, '--freq', quiet = TRUE)

  plink_frq <- 'plink.frq'

  expect_true(file.exists(plink_frq))
  df <- read_plink_freq(plink_frq)
  unlink(plink_frq)
  expect_equal(nrow(df), bed_nb_snps(bo))

  ### snp subset
  bed_plink_cmd(bo, snp_idx = 10:6, '--freq', quiet = TRUE)
  df2 <- read_plink_freq('plink.frq')
  bim_df <- bed_bim_df(bo, subset = FALSE)
  # N.B: the snp ordering is lost
  expect_identical(df2$SNP, bim_df[6:10, 'SNP'])

  idx <- c(2, 10:16)
  bo2 <- bed_subset_snps_by_idx(bo, idx)

  bed_plink_cmd(bo2, '--freq', quiet = TRUE)
  df2 <- read_plink_freq('plink.frq')
  expect_identical(df2$SNP, bim_df[idx, 'SNP'])

  ### sample subset
  sample_idx <- c(3, 11:20, 50:80)
  bed_plink_cmd(bo, sample_idx = sample_idx, '--freq', quiet = TRUE)
  df2 <- read_plink_freq('plink.frq')
  expect_equal(max(df2$NCHROBS), length(sample_idx) * 2)

  bo2 <- bed_subset_samples_by_idx(bo, sample_idx)
  bed_plink_cmd(bo2, '--freq', quiet = TRUE)
  df3 <- read_plink_freq('plink.frq')
  expect_identical(df3, df2)


  ### both subsets
  bo2 <- bed_subset_snps_by_idx(bo, idx)
  bo3 <- bed_subset_samples_by_idx(bo2, sample_idx)
  bed_plink_cmd(bo3, '--freq', quiet = TRUE)
  df4 <- read_plink_freq('plink.frq')

  expect_identical(df4$SNP, bim_df[idx, 'SNP'])
  expect_equal(max(df4$NCHROBS), length(sample_idx) * 2)

  ### missing_phenotype
  expect_error(bed_plink_cmd(bo, '--freq', missing_phenotype = 1:2, quiet = TRUE),
    'bad missing_phenotype')
  expect_error(bed_plink_cmd(bo, '--freq', missing_phenotype = 1, quiet = TRUE),
    'bad missing_phenotype')
  expect_error(bed_plink_cmd(bo, '--freq', missing_phenotype = NA_integer_,
      quiet = TRUE), 'bad missing_phenotype')

  ### alleles order
  bed_plink_cmd(bo, '--freq', keep_allele_order = TRUE, quiet = TRUE)
  df1 <- read_plink_freq('plink.frq')

  expect_error(bed_plink_cmd(bo, '--freq', a2_alleles = 'A', quiet = TRUE),
    'bad a2_alleles length')

  a2 <- bed_allele_higher(bo)
  bed_plink_cmd(bo, '--freq', a2_alleles = a2, keep_allele_order = TRUE,
    quiet = TRUE)
  df2 <- read_plink_freq('plink.frq')
  expect_false(all(df1$MAF == df2$MAF))

  inv <- which(df1$A1 != df2$A1)
  df2$MAF[inv] <- 1 - df2$MAF[inv]
  df2$A1[inv] <- df1$A1[inv]
  df2$A2[inv] <- df1$A2[inv]
  expect_equal(df2, df1, tolerance = 1e-4)

  df2 <- read_plink_freq('plink.frq') # reread

  bed_plink_cmd(bo, '--freq', lexicographic_allele_order = TRUE, quiet = TRUE)
  df3 <- read_plink_freq('plink.frq') # reread
  expect_identical(df3 ,df2)

}
test_that('bed_plink_cmd', .bed_plink_cmd())



.plink_cmd <- function() {
  plink_cmd <- plinker:::plink_cmd

  # edge cases
  suppressWarnings(
    expect_error(plink_cmd('toto', stderr = TRUE, stdout = TRUE),
      'ERROR running plink'))

  suppressWarnings(
    expect_error(plink_cmd('toto', command = '/not/found/foo',
         stderr = TRUE, stdout = TRUE),
      'error in running command'))

  text <- plink_cmd('--help', stdout = TRUE)
  expect_gt(length(text), 100)
}
test_that('plink_cmd', .plink_cmd())




.plink_version <- function() {
  ver <- plink_version()
  expect_match(ver, "PLINK v")
}
test_that('plink_version', .plink_version())

