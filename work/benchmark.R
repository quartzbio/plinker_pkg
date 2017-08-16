library(devtools)
load_all('plinker')


###############################
# BEDMatrix
############################
system.time(bo <- bed_open('1kg_phase1_chr1'))
library(BEDMatrix)
bmat <- BEDMatrix(bo$bed, n = bed_nb_samples(bo, subset = FALSE),
  p = bed_nb_snps(bo, subset = FALSE))

bmat
# BEDMatrix: 1092 x 3007196 [/home/docker/workspace/plinker_pkg/1kg_phase1_chr1.bed]
str(bmat)
# BEDMatrix: 1092 x 3007196 [/home/docker/workspace/plinker_pkg/1kg_phase1_chr1.bed]

dim(bmat)
# [1]    1092 3007196
length(bmat)
# [1] 3283858032
dimnames(bmat)
# [[1]]
# NULL
#
# [[2]]
# NULL
dimnames(bmat)

bmat[1:3, 1:5]
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    0    0    0    0    0
# [2,]    0    1    1    1    1
# [3,]    0    0    0    0    0

bmat[c(2, 3), 1:5]
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    0    1    1    1    1
# [2,]    0    0    0    0    0
bmat[c(3, 2), 1:5]
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    0    0    0    0    0
# [2,]    0    1    1    1    1

bmat[1, 1:5]
# [1] 0 0 0 0 0
bmat[1, 1:5, drop = FALSE]
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    0    0    0    0    0


path <- system.file("extdata", "example.bed", package = "BEDMatrix")
m <- BEDMatrix(path)


m[1:3, c("snp0_A", "snp1_C", "snp2_G")]
#           snp0_A snp1_C snp2_G
# per0_per0      0      1      1
# per1_per1      1      1      1
# per2_per2      1      0      0

m[1:3, c("snp0_", "snp1_C", "snp2_G")]
# Error in convertIndex(x, j, "j") : subscript out of bounds



# p

2L - bmat[1:3, 1:5]




# toto.fam has ~ 100000 lines/samples

system.time(df <- read_fam('toto.fam'))
#    user  system elapsed
#   0.069   0.000   0.067

#    user  system elapsed
#   0.071   0.000   0.070
dim(df)
# [1] 100036      6


system.time(df <- read_fam('1kg_phase1_all.fam'))

system.time(df <- read_bim('1kg_phase1_all.bim'))
#    user  system elapsed
#  70.496   1.331  71.949

# user  system elapsed
# 174.905   1.463 176.706

# test only read ID column
path <- '1kg_phase1_all.bim'
system.time(
df <- as.data.frame(readr::read_delim(path, '\t',
    col_names = c('CHR', 'SNPID', 'MORGANS', 'POS', 'A1', 'A2'),
    col_types = readr::cols_only(SNPID = 'c'),
    progress = FALSE
  ))
)
#    user  system elapsed
# 105.671   0.286 106.046
#user  system elapsed
#70.496   1.331  71.949

# test with skip
system.time(
  df <- as.data.frame(readr::read_delim(path, '\t',
          col_names = c('CHR', 'SNPID', 'MORGANS', 'POS', 'A1', 'A2'),
          col_types = readr::cols_only(SNPID = 'c'),
          progress = FALSE,
          skip = 2e7
      ))
)
#    user  system elapsed
#  14.122   0.195  14.339



### using data.table::fread: SNPIDS
path <- '1kg_phase1_all.bim'
system.time(
df <- data.table::fread(path, sep = '\t', header = FALSE, verbose = FALSE,
  select = 2, quote = '', strip.white = FALSE, data.table = FALSE,
  colClasses = rep('character', 6))
)

#  15.572   0.647  16.246

#    user  system elapsed
#  15.641   0.415  16.078

path <- '1kg_phase1_chr1.bim'
system.time(
  df <- data.table::fread(path, sep = '\t', header = FALSE, verbose = FALSE,
      select = 2, quote = '', strip.white = FALSE, data.table = FALSE,
      colClasses = rep('character', 6))
)
#    user  system elapsed
#   1.139   0.032   1.173



### using data.table::fread: all
path <- '1kg_phase1_all.bim'
system.time(
  df <- data.table::fread(path, sep = '\t', header = FALSE, verbose = FALSE,
      quote = '', strip.white = FALSE, data.table = FALSE)
)

#    user  system elapsed
#  23.254   0.692  23.977

#  33.410   0.973  34.460


system.time(df <- read_bim('1kg_phase1_chr1.bim'))
#    user  system elapsed
#   9.256   0.110   9.380

system.time(lines <- readLines('1kg_phase1_chr1.bim'))
#    user  system elapsed
#   8.901   0.000   8.849

system.time(lines <- data.table::fread('1kg_phase1_chr1.bim'))
#    user  system elapsed
#   1.622   0.032   1.656



readLines('1kg_phase1_chr1.bim', 2)
# [1] "1\trs58108140\t0\t10583\tA\tG"  "1\trs189107123\t0\t10611\tG\tC"
con <- file('1kg_phase1_chr1.bim')
readLines(con, 1)
# [1] "1\trs58108140\t0\t10583\tA\tG"
readLines(con, 1)
# [1] "1\trs58108140\t0\t10583\tA\tG"

readLines




system.time(bo <- bed_open('1kg_phase1_chr1'))

system.time(bim_df <- bed_bim_df(bo))
#    user  system elapsed
#       0       0       0

#    user  system elapsed
#   2.599   0.074   2.678


system.time(bo <- bed_open('1kg_phase1_all'))

bed_nb_snps(bo)
# [1] 39728178



### split_sorted_ints_by_blocks
# worst case
idx <- seq.int(1, 39728178, by = 2)

system.time(blks <- plinker:::split_sorted_ints_by_blocks(idx))
#    user  system elapsed
#   1.607   0.304   1.915

length(idx)/nrow(blks)
# [1] 1


# best case
idx <- c(1:100000, 1e7, 3e7:39728178)
system.time(blks <- plinker:::split_sorted_ints_by_blocks(idx))
#    user  system elapsed
#   0.356   0.049   0.406

length(idx)/nrow(blks)
# [1] 1

#####################################################

system.time(bo <- bed_open('1kg_phase1_chr1'))
bo2 <- bed_subset_snps_by_idx(bo, 1001:2000)


bof <- bed_filter_snps_by_missing_rate(bo, 0.000001)


system.time(counts <- bed_plink_freq_count(bo2))
#    user  system elapsed
#   1.014   0.164   1.161

system.time(counts <- bed_plink_freq_count(bo))
#    user  system elapsed
#  28.739   0.625  20.183

bads <- with(counts, which(C1 > C2))
head(counts[bads, ])
#    CHR         SNP A1 A2   C1  C2 G0
# 7    1 rs140337953  T  G 1570 596  0
# 18   1 rs150021059  G  T 1929 237  0
# 25   1   rs3091274  A  C 1941 225  0
# 40   1 rs189727433  C  A 1887 279  0
# 44   1  rs74970982  G  A 2067  99  0
# 67   1  rs75062661  G  A 1414 752  0

bim_df <- bed_bim_df(bo, subset = FALSE)
head(bim_df[bads, ])
#    CHR       SNPID MORGANS   POS A1 A2
# 7    1 rs140337953       0 30923  T  G
# 18   1 rs150021059       0 52238  G  T
# 25   1   rs3091274       0 55164  A  C
# 40   1 rs189727433       0 57952  C  A
# 44   1  rs74970982       0 61442  G  A
# 67   1  rs75062661       0 69511  G  A

library(dplyr)
nb <- unique(with(counts, C1 + C2 + 2*G0) / 2)
nb
# [1] 1092

bed_nb_samples(bo)
# [1] 1092




