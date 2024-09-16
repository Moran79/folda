test_that("1D plot with >2 classes", {
  skip_on_cran()
  dat <- iris[, c(5,1)]
  fit <- folda(dat[, -1, drop = FALSE], dat[,1])
  plot(fit, dat, dat[,1]) # check if this plot works
  expect_equal(2 * 2, 4)
})

test_that("1D plot: class with only 1 obs", {
  skip_on_cran()
  dat <- iris[c(1, 51:60, 101:110), c(5,1)]
  fit <- folda(dat[, -1, drop = FALSE], dat[,1])
  plot(fit, dat, dat[,1]) # check if this plot works
  expect_equal(2 * 2, 4)
})

test_that("1D plot: dominant class", {
  skip_on_cran()
  dat <- iris[c(1:52, 101:102), c(5,2)]
  fit <- folda(dat[, -1, drop = FALSE], dat[,1])
  plot(fit, dat, dat[,1]) # check if this plot works
  expect_equal(2 * 2, 4)
})
