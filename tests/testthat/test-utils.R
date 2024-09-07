test_that("Downsampling and Balanced Group", {
  skip_on_cran()
  response <- factor(c("A", "A", NA, "A", "B", "B", "C", "C", "C", NA), levels = LETTERS[1:5])
  set.seed(443)
  result <- getDownSampleInd(response, downSampling = TRUE)
  result2 <- getDownSampleInd(response, downSampling = TRUE, kSample = 1)
  expected_result <- c(4, 2, 5, 6, 9, 8)
  expected_result2 <- c(1, 5, 7)
  expect_equal(result, expected_result)
  expect_equal(result2, expected_result2)
})


test_that("Missing Value Imputation", {
  skip_on_cran()
  dat <- data.frame(X1 = rep(NA,5),
                    X2 = factor(rep(NA,5), levels = LETTERS[1:3]),
                    X3 = 1:5,
                    X4 = LETTERS[1:5],
                    X5 = c(NA, 2, 3, 10, NA),
                    X6 = factor(c("A", NA, NA, "B", "B"), levels = LETTERS[1:3]))
  result1 <- missingFix(dat, missingMethod = c("medianFlag", "mode"))
  result2 <- missingFix(dat, missingMethod = c("mean", "newLevel"))
  expected_result1 <- data.frame(X3 = 1:5,
                                 X4 = factor(LETTERS[1:5]),
                                 X5 = c(3, 2, 3, 10, 3),
                                 X6 = factor(LETTERS[c(1,2,2,2,2)]),
                                 X5_FLAG = c(1, 0, 0, 0, 1))
  expected_result1ref <- data.frame(X3 = 3,
                                    X4 = factor("A", levels = LETTERS[1:5]),
                                    X5 = 3,
                                    X6 = factor("B", levels = LETTERS[1:2]),
                                    X5_FLAG = 1)
  expected_result2 <- data.frame(X3 = 1:5,
                                 X4 = factor(LETTERS[1:5]),
                                 X5 = c(5, 2, 3, 10, 5),
                                 X6 = factor(c("A", "new0_0Level", "new0_0Level", "B", "B")))
  expected_result2ref <- data.frame(X3 = 3,
                                    X4 = factor("A", levels = LETTERS[1:5]),
                                    X5 = 5,
                                    X6 = factor("new0_0Level", levels = c("A", "B", "new0_0Level")))
  expect_equal(result1$data, expected_result1)
  expect_equal(result1$ref, expected_result1ref)
  expect_equal(result2$data, expected_result2)
  expect_equal(result2$ref, expected_result2ref)
})


test_that("getDataInShape", {
  skip_on_cran()
  dat <- data.frame(X1_FLAG = c(0, 0, 0),
                    X1 = factor(c(NA, "C", "B"), levels = LETTERS[2:3]),
                    X2_FLAG = c(NA, 0, 1),
                    X4 = 1:3)
  misRef <- data.frame(X1_FLAG = 1,
                       X1 = factor("A", levels = LETTERS[1:2]),
                       X2_FLAG = 1,
                       X3 = factor("B", levels = LETTERS[1:2]),
                       X3_FLAG = 1,
                       X4 = 1,
                       X4_FLAG = 1)
  result <- getDataInShape(dat, misRef)
  expected_result <- data.frame(X1_FLAG = c(1, 1, 0),
                                X1 = factor(c("A", "A", "B"), levels = LETTERS[1:2]),
                                X2_FLAG = c(1, 0, 1),
                                X3 = factor(c("B", "B", "B"), levels = LETTERS[1:2]),
                                X3_FLAG = c(1, 1, 1),
                                X4 = 1:3,
                                X4_FLAG = c(0, 0, 0))
  expect_equal(result, expected_result)
})


