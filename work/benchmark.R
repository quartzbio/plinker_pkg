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