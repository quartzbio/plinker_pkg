context('distances')

.ibs_state <- function() {
  ibs_state <- plinker:::ibs_state

  expect_equal(ibs_state(0, 0), 2)
  expect_equal(ibs_state(2, 2), 2)

  expect_equal(ibs_state(1, 1), 2)

  expect_equal(ibs_state(0, 1), 1)
  expect_equal(ibs_state(1, 0), 1)

  expect_equal(ibs_state(0, 2), 0)
  expect_equal(ibs_state(2, 0), 0)

  expect_equal(ibs_state(0, NA), NA_integer_)

  expect_equal(ibs_state(
      c(0,1,2,NA,2,1,0,1,2,0),
      c(0,1,2,1 ,0,1,2,0,1,2)),
      c(2,2,2,NA_integer_,0,2,0,1,1,0))

}
test_that('ibs_state', .ibs_state())

.ibs <- function() {
  ibs <- plinker:::ibs

  expect_equal(ibs(0), 0)
  expect_equal(ibs(1), 0.5)
  expect_equal(ibs(2), 1)

  expect_equal(ibs(rep(2L, 10)), 1)
  expect_equal(ibs(rep(1L, 10)), 0.5)
  expect_equal(ibs(rep(0L, 10)), 0)

  expect_equal(ibs(c(0,1,2,NA)), 3/6)
  expect_equal(ibs(NA_integer_), NA_real_)

}
test_that('ibs', .ibs())