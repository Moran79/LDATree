test_that("folda: work on tibble", {
  skip_on_cran()
  dat <- ggplot2::diamonds[1:100,]
  fit <- Treee(dat[, -2], response = dat[[2]], verbose = FALSE)
  result <- predict(fit, dat)
  expect_equal(result[1:4], c("Very Good", "Ideal", "Ideal", "Premium"))
})
