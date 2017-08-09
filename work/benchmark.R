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



