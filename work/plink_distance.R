library(devtools)
load_all('plinker')

bo <- bed_open('1kg_phase1_chr1')


system.time(dist <- bed_plink_distance(bo))
#    user  system elapsed
# 216.627   5.313  58.190

dim(dist)
# [1] 1092 1092
sum(dist)
# [1] 1122610

bo2 <- bed_subset(bo, snp_idx = 1:1000000)
#

system.time(dist <- bed_plink_distance(bo2, threads = 1))
#    user  system elapsed
#  24.879   0.223  25.069
system.time(dist <- bed_plink_distance(bo2, threads = 2))
#    user  system elapsed
#  29.888   0.384  18.012

system.time(dist <- bed_plink_distance(bo2, threads = 3))

system.time(dist <- bed_plink_distance(bo2, threads = 4))
#    user  system elapsed
#  37.143   0.586  15.490


system.time(dist <- bed_plink_distance(bo2, threads = 8))
#    user  system elapsed
#  56.565   0.907  17.336



system.time(dist <- bed_plink_distance(bo2))
#    user  system elapsed
#  75.061   2.010  19.730





####### NOT MULTITHREADED !! #######

system.time(df1 <- bed_plink_fisher(bo2, phenotype = pheno, threads = 1))


