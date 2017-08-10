library(devtools)

check_man('plinker')

test('plinker')
test('plinker', 'bed')
test('plinker', 'bim')
test('plinker', 'fam')
test('plinker', 'filters')
test('plinker', 'plink')
test('plinker', 'sample_data')
test('plinker', 'snp_subset')
test('plinker', 'sample_subset')
test('plinker', 'utils')