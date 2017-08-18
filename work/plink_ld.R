library(devtools)
load_all('plinker')

bo <- bed_open('1kg_phase1_chr1')
system.time(dfs <- bed_plink_ld(bo, threads = 1, window_size = 2))
#    user  system elapsed
#  30.016   1.154  31.207

system.time(ld2 <- bed_plink_ld(bo, threads = 8, window_size = 2))
#    user  system elapsed
#  30.370   1.309  30.568


### not multithreaded

system.time(ld1 <- bed_plink_ld(bo, threads = 1, window_size = 2, min_r2 = 0.9))
#    user  system elapsed
#  26.736   0.898  27.665
nrow(ld1)
# [1] 63472

system.time(ld1 <- bed_plink_ld(bo, threads = 8, window_size = 2, min_r2 = 0.9))
#    user  system elapsed
#  27.162   0.980  26.891

### allele order impact
system.time(ld1 <- bed_plink_ld(bo, window_size = 2, keep_allele_order = FALSE))
#    user  system elapsed
#  31.742   1.275  30.599

system.time(ld2 <- bed_plink_ld(bo, window_size = 2, keep_allele_order = TRUE))
#    user  system elapsed
#  31.308   1.293  30.413
