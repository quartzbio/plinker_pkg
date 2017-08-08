library(devtools)
load_all('plinker')

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




system.time(df <- read_bim('1kg_phase1_chr1.bim'))
#    user  system elapsed
#   9.256   0.110   9.380


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
#   2.381   0.347   2.731
nrow(blks)
# [1] 19864089
nrow(blks)/length(idx)
# [1] 1


