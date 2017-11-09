library(devtools)
load_all('plinker')



bo <- bed_open('1kg_phase1_chr1')
system.time(dfs <- bed_plink_missing(bo))
#    user  system elapsed
#   7.155   0.615   9.129
head(dfs$snp)
#   CHR         SNP N_MISS N_GENO F_MISS
# 1   1  rs58108140      0   1092      0
# 2   1 rs189107123      0   1092      0
# 3   1 rs180734498      0   1092      0
# 4   1 rs144762171      0   1092      0
# 5   1 rs201747181      0   1092      0
# 6   1 rs151276478      0   1092      0

head(dfs$sample)
#   FID     IID MISS_PHENO N_MISS  N_GENO F_MISS
# 1   0 HG00096          Y      0 3007196      0
# 2   0 HG00097          Y      0 3007196      0
# 3   0 HG00099          Y      0 3007196      0
# 4   0 HG00100          Y      0 3007196      0
# 5   0 HG00101          Y      0 3007196      0
# 6   0 HG00102          Y      0 3007196      0

