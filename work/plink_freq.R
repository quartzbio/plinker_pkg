library(devtools)
load_all('plinker')


bo <- bed_open('1kg_phase1_chr1')
system.time(counts <- bed_plink_freq_count(bo))
#    user  system elapsed
#   3.546   0.428   3.956
head(counts)
#   CHR         SNP A1 A2  C1   C2 G0
# 1   1  rs58108140  A  G 314 1870  0
# 2   1 rs189107123  G  C  41 2143  0
# 3   1 rs180734498  T  C 249 1935  s
# 4   1 rs144762171  C  G  59 2125  0
# 5   1 rs201747181  T TC  35 2149  0
# 6   1 rs151276478  C  T  45 2139  0


all(counts$A1 == bed_allele1(bo))
# [1] TRUE


system.time(freqs <- bed_plink_freqx(bo))
#    user  system elapsed
#   4.017   0.382   4.388
head(freqs)
#   CHR         SNP A1 A2 C(HOM A1) C(HET) C(HOM A2) C(HAP A1) C(HAP A2)
# 1   1  rs58108140  A  G         5    304       783         0         0
# 2   1 rs189107123  G  C         0     41      1051         0         0
# 3   1 rs180734498  T  C         6    237       849         0         0
# 4   1 rs144762171  C  G         0     59      1033         0         0
# 5   1 rs201747181  T TC         1     33      1058         0         0
# 6   1 rs151276478  C  T         0     45      1047         0         0
#   C(MISSING)
# 1          0
# 2          0
# 3          0
# 4          0
# 5          0
# 6          0
all(freqs$A1 == bed_allele1(bo))
# [1] TRUE

