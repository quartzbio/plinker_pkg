library(devtools)
load_all('plinker')

bo <- bed_open('1kg_phase1_chr1')
nbind <- bed_nb_samples(bo)
nbind
# [1] 1092
pheno <- c(rep(1L, nbind/2), rep(2L, nbind/2))


system.time(df1 <- bed_plink_fisher(bo, phenotype = pheno, keep_allele_order = FALSE))
#    user  system elapsed
#  37.788   2.112  40.271

system.time(df2 <- bed_plink_fisher(bo, phenotype = pheno, keep_allele_order = TRUE))
#    user  system elapsed
#  36.220   2.187  39.043

identical(df2, df1)
# [1] FALSE

bim <- bed_bim_df(bo)
bim[7, ]
#   CHR         SNP MORGANS   POS A1 A2
# 7   1 rs140337953       0 30923  T  G


df1[31:35, ]
#    CHR         SNP A1 A2    TEST        AFF      UNAFF         P
# 31   1 rs140337953  G  T    GENO 91/178/277 55/130/361 1.096e-06
# 32   1 rs140337953  G  T   TREND    360/732    240/852 4.184e-07
# 33   1 rs140337953  G  T ALLELIC    360/732    240/852 1.074e-08
# 34   1 rs140337953  G  T     DOM    269/277    185/361 3.269e-07
# 35   1 rs140337953  G  T     REC     91/455     55/491 1.784e-03

df2[31:35, ]
#    CHR         SNP A1 A2    TEST        AFF      UNAFF         P
# 31   1 rs140337953  T  G    GENO 277/178/91 361/130/55 1.096e-06
# 32   1 rs140337953  T  G   TREND    732/360    852/240 4.184e-07
# 33   1 rs140337953  T  G ALLELIC    732/360    852/240 1.074e-08
# 34   1 rs140337953  T  G     DOM     455/91     491/55 1.784e-03
# 35   1 rs140337953  T  G     REC    277/269    361/185 3.269e-07






####### NOT MULTITHREADED !! #######

system.time(df1 <- bed_plink_fisher(bo, phenotype = pheno, threads = 1))
#    user  system elapsed
#  44.457   2.041  46.948


system.time(df2 <- bed_plink_fisher(bo, phenotype = pheno, threads = 8))
#    user  system elapsed
#  39.103   2.496  41.987

identical(df2, df1)
# [1] TRUE

## recompute pvalues
library(dplyr)
df <- filter(df1, TEST == 'ALLELIC') %>% head

aff <- strsplit(df$AFF, '/')
aff <- lapply(aff, as.integer)
unaff <- strsplit(df$UNAFF, '/')
unaff <- lapply(unaff, as.integer)


cts <- mapply(rbind, aff, unaff, SIMPLIFY = FALSE)

res <- lapply(cts, fisher.test)
pvs <- unlist(lapply(res, getElement, 'p.value'))

df$P
# [1] 4.941e-04 1.138e-01 1.496e-05 1.125e-01 3.066e-01 7.636e-01
pvs
# [1] 4.941473e-04 1.138327e-01 1.495712e-05 1.125115e-01 3.065787e-01
# [6] 7.635989e-01

abs(df$P - pvs) / pmax(df$P, pvs)
# [1] 9.573981e-05 2.874499e-04 1.924838e-04 1.018080e-04 6.937174e-05
# [6] 1.402181e-06
# ==> seems to match

