library(devtools)

check_man(as_pkg('plinker'))

test('plinker')
test('plinker', 'bed')
test('plinker', 'sample_data')

