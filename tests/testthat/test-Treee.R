test_that("folda: work on tibble", {
  skip_on_cran()
  dat <- ggplot2::diamonds[1:100,]
  fit <- Treee(dat[, -2], response = dat[[2]], verbose = FALSE)
  result <- predict(fit, dat)
  expect_equal(result[1:4], c("Ideal", "Premium", "Premium", "Very Good"))
})

test_that("folda: all columns are constant", {
  skip_on_cran()
  dat <- iris[, c(5,1:4)]; dat[, 2:5] <- dat[rep(1,150), 2:5]
  fit <- Treee(dat[, -1], dat[, 1], verbose = FALSE)
  expect_equal(fit[[1]]$nodeModel, "mode")
})
