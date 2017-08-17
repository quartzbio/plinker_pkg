library(devtools)

check_man('plinker')

test('plinker')
test('plinker', 'bed$')
test('plinker', 'bedmatrix')
test('plinker', 'bim')
test('plinker', 'fam')
test('plinker', 'filters')
test('plinker', 'plink')
test('plinker', 'sample_data')
test('plinker', 'snp_subset')
test('plinker', 'sample_subset')
test('plinker', 'utils')

help2 <- function(..., help_type = 'html') {
  doc <- utils::help(..., help_type = help_type)
  utils:::print.help_files_with_topic(doc)
}
