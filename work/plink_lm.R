library(devtools)
load_all('plinker')

mk_pheno <- function(bo) {
  set.seed(123)
  rnorm(bed_nb_samples(bo))
}

bo <- bed_open('1kg_phase1_chr1')


freqx <- bed_plink_freqx(bo)
non_rare <- which(freqx$`C(HOM A1)` > 100)

bo2 <- bed_subset(bo, snp_idx = non_rare)


bo1 <- bed_subset(bo, snp_idx= 1:2, sample_idx = 1:10)
pheno1 <- mk_pheno(bo1)
fam1 <- bed_fam_df(bo1)


bed_plink_lm(bo1, pheno1)

covars1 <- fam1[, 1:2]
covars1$TOTO <- c(rep(1, 5), c(rep(2, 5)))
covars1$TITI <- 1:nrow(covars1)
res1 <- bed_plink_lm(bo1, pheno1, covars = covars1)
#   CHR         SNP    BP A1 TEST NMISS    BETA    STAT      P
# 1   1  rs58108140 10583  A  ADD    10  0.4969  0.6584 0.5347
# 2   1  rs58108140 10583  A TOTO    10  0.6100  0.4423 0.6738
# 3   1  rs58108140 10583  A TITI    10 -0.1895 -0.7939 0.4575
# 4   1 rs189107123 10611  G  ADD    10 -0.8176 -0.6638 0.5315
# 5   1 rs189107123 10611  G TOTO    10  0.7502  0.5470 0.6041
# 6   1 rs189107123 10611  G TITI    10 -0.2303 -0.9350 0.3859

res2 <- bed_R_lm(bo1, pheno1, covars = covars1)
#   CHR         SNP    BP A1 TEST NMISS       BETA       STAT         P
# 1   1  rs58108140 10583  A  ADD    10  0.4968941  0.6584306 0.5346955
# 2   1 rs189107123 10611  G TOTO    10  0.6099845  0.4422743 0.6737872
# 3   1  rs58108140 10583  A TITI    10 -0.1894505 -0.7938564 0.4575172
# 4   1 rs189107123 10611  G  ADD    10 -0.8175977 -0.6638060 0.5314785
# 5   1  rs58108140 10583  A TOTO    10  0.7502432  0.5470069 0.6041070
# 6   1 rs189107123 10611  G TITI    10 -0.2303304 -0.9350241 0.3858627
all.equal(res2, res, tolerance = 1e-4)

covars1 <- fam1[, 1:2]
covars1$TOTO <-  LETTERS[1:nrow(covars1)]
bed_plink_lm(bo1, pheno1, covars = covars1)
# Error in read_plink_output("plink.assoc.linear") :
#   plink output file does not exist:plink.assoc.linear


bed_plink_lm(bo1, pheno1)
#   CHR         SNP    BP A1 TEST NMISS    BETA    STAT      P
# 1   1  rs58108140 10583  A  ADD    10  0.4166  0.6105 0.5585
# 2   1 rs189107123 10611  G  ADD    10 -0.3387 -0.3196 0.7574



df1 <- bed_R_lm(bo1, pheno1)
df1
#   CHR         SNP MORGANS   POS A1 A2        BETA      STAT      PVAL
# 1   1  rs58108140       0 10583  A  G -0.06726087 1.0605111 0.3033271
# 2   1 rs189107123       0 10611  G  C -0.11446071 0.5182898 0.4717264

pheno2
df2 <- bed_plink_lm(bo1, pheno1)
df2
#   CHR         SNP    BP A1 TEST NMISS     BETA    STAT      P
# 1   1  rs58108140 10583  A  ADD  1092 -0.06726 -1.0300 0.3033
# 2   1 rs189107123 10611  G  ADD  1092 -0.11450 -0.7199 0.4717


good_samples <- which(!is.na(pheno))
bo3 <- bed_subset(bo1, sample_idx = good_samples)

pheno2 <- pheno1
pheno2[1] <- -9L
df2 <- bed_plink_lm(bo1, pheno2)
df2
#   CHR         SNP    BP A1 TEST NMISS   BETA    STAT      P
# 1   1  rs58108140 10583  A  ADD  1091 -0.068 -1.0410 0.2982
# 2   1 rs189107123 10611  G  ADD  1091 -0.115 -0.7232 0.4697



bb <- bed_subset(bo, snp_idx = c(57, 154:100))
bbb <- bed_subset(bb, snp_idx = 2)
bed_snp_idx(bbb)
# [1] 100

bed_bim_df(bbb)
#     CHR         SNP MORGANS   POS A1 A2
# 100   1 rs183605470       0 84139  T  A
head(bed_bim_df(bbb, subset = FALSE))
#   CHR         SNP MORGANS   POS A1 A2
# 1   1  rs58108140       0 10583  A  G
# 2   1 rs189107123       0 10611  G  C
# 3   1 rs180734498       0 13302  T  C
# 4   1 rs144762171       0 13327  C  G
# 5   1 rs201747181       0 13957  T TC
# 6   1 rs151276478       0 13980  C  T

bed_subset(bo, snp_IDs = 'rs180734498')




out <- bed_plink_lm(bo2, mk_pheno(bo2))
head(out)
#   CHR         SNP    BP A1 TEST NMISS    BETA   STAT       P
# 1   1 rs140337953 30923  G  ADD  1092 0.06514 1.5480 0.12190
# 2   1 rs150021059 52238  T  ADD  1092 0.06693 1.1160 0.26480
# 3   1   rs3091274 55164  C  ADD  1092 0.04329 0.7055 0.48060
# 4   1 rs189727433 57952  A  ADD  1092 0.11590 2.1100 0.03507
# 5   1  rs74970982 61442  A  ADD  1092 0.18870 1.9190 0.05524
# 6   1 rs201888535 63735  C  ADD  1092 0.01260 0.3154 0.75250

out2 <- bed_plink_lm(bo2, mk_pheno(bo2), model = 'genotypic')
head(out2)
#   CHR         SNP    BP A1     TEST NMISS     BETA    STAT       P
# 1   1 rs140337953 30923  G      ADD  1092  0.08042  1.7570 0.07924
# 2   1 rs140337953 30923  G   DOMDEV  1092 -0.06191 -0.8481 0.39660
# 3   1 rs140337953 30923  G GENO_2DF  1092       NA  3.1150 0.21060
# 4   1 rs150021059 52238  T      ADD  1092  0.14630  1.8990 0.05782
# 5   1 rs150021059 52238  T   DOMDEV  1092 -0.18240 -1.6410 0.10120
# 6   1 rs150021059 52238  T GENO_2DF  1092       NA  3.9380 0.13960

out3 <- bed_plink_lm(bo2, mk_pheno(bo2), model = 'hethom')
head(out3)
#   CHR         SNP    BP A1     TEST NMISS     BETA    STAT       P
# 1   1 rs140337953 30923  G      HOM  1092  0.16080  1.7570 0.07924
# 2   1 rs140337953 30923  G      HET  1092  0.01851  0.2673 0.78930
# 3   1 rs140337953 30923  G GENO_2DF  1092       NA  3.1150 0.21060
# 4   1 rs150021059 52238  T      HOM  1092  0.29260  1.8990 0.05782
# 5   1 rs150021059 52238  T      HET  1092 -0.03607 -0.4155 0.67780
# 6   1 rs150021059 52238  T GENO_2DF  1092       NA  3.9380 0.13960

out4 <- bed_plink_lm(bo2, mk_pheno(bo2), model = 'dominant')
head(out4)
#   CHR         SNP    BP A1 TEST NMISS    BETA   STAT       P
# 1   1 rs140337953 30923  G  DOM  1092 0.06428 1.0490 0.29460
# 2   1 rs150021059 52238  T  DOM  1092 0.03659 0.4673 0.64030
# 3   1   rs3091274 55164  C  DOM  1092 0.01225 0.1526 0.87870
# 4   1 rs189727433 57952  A  DOM  1092 0.11750 1.5660 0.11770
# 5   1  rs74970982 61442  A  DOM  1092 0.20000 1.8580 0.06341
# 6   1 rs201888535 63735  C  DOM  1092 0.02341 0.3841 0.70100

out5 <- bed_plink_lm(bo2, mk_pheno(bo2), model = 'recessive')
head(out5)
#   CHR         SNP    BP A1 TEST NMISS     BETA   STAT       P
# 1   1 rs140337953 30923  G  REC  1092 0.154800 1.7450 0.08120
# 2   1 rs150021059 52238  T  REC  1092 0.297900 1.9410 0.05249
# 3   1   rs3091274 55164  C  REC  1092 0.237000 1.5090 0.13150
# 4   1 rs189727433 57952  A  REC  1092 0.306600 2.3170 0.02068
# 5   1  rs74970982 61442  A  REC  1092 0.373700 0.9141 0.36090
# 6   1 rs201888535 63735  C  REC  1092 0.009298 0.1219 0.90300

