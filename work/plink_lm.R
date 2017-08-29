library(devtools)
load_all('plinker')

bo <- bed_open('1kg_phase1_chr1')

miss <- bed_plink_missing(bo)

freqx <- bed_plink_freqx(bo)
non_rare <- which(freqx$`C(HOM A1)` > 100)

bo2 <- bed_subset(bo, snp_idx = non_rare)


mk_pheno <- function(bo) {
  set.seed(123)
  rnorm(bed_nb_samples(bo))
}

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

