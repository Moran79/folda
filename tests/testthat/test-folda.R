test_that("Forward Selection on Dummy", {
  skip_on_cran()
  set.seed(443)
  dat <- data.frame(Class = as.factor(sample(LETTERS[1:4],1000,replace = TRUE)))
  dat <- cbind(dat, stats::model.matrix(~.-1, data = dat), Z = matrix(rnorm(4000),1000,4))
  fitPillai <- folda(dat[,-1], response = dat[,1])
  fitWilks <- folda(dat[,-1], response = dat[,1], testStat = "Wilks")

  expect_equal(length(fitPillai$forwardInfo$var), 3)
  expect_equal(length(fitWilks$forwardInfo$var), 1)
})

test_that("folda: Prior and Misclassification Cost", {
  skip_on_cran()
  fitPillai <- folda(iris[, -5],response = iris[, 5],prior = c(0,1,2),
                     misClassCost = matrix(c(0, 1, 0,
                                             1, 0, 100,
                                             1, 1, 0), 3, 3, byrow = TRUE))
  result <- as.vector(table(stats::predict(fitPillai, iris), dnn = NULL))
  expect_equal(result, c(61L, 89L))
})

test_that("folda: Correction and alpha", {
  skip_on_cran()
  fitPillai <- folda(iris[, -5], response = iris[, 5], correction = FALSE, alpha = 0.8)
  result <- round(fitPillai$forwardInfo$threshold, 5)
  expect_equal(result, c(0.00303, 0.00053, 0.00027, 0.00019))
})
