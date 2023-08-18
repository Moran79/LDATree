test_that("works on iris data", {
  expect_equal(predict(Treee(Species~., data = iris),iris)[c(1,51,101)], c("setosa", "versicolor", "virginica"))
})
