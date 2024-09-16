test_that("updatePriorAndMisClassCost works", {
  skip_on_cran()

  # Simple case
  test1 <- updatePriorAndMisClassCost(prior = NULL,
                                      misClassCost = NULL,
                                      response = factor(LETTERS[c(1,1,1,2,2,3)]),
                                      insideNode = FALSE)
  expect1 <- list(prior = c(A = 1, B = 1, C = 1), misClassCost = structure(
    c(0, 1, 1, 1, 0, 1, 1, 1, 0), dim = c(3L, 3L), dimnames = list(c("A", "B", "C"), c("A", "B", "C"))))
  expect_equal(test1, expect1)

  # prior != obs
  test2 <- updatePriorAndMisClassCost(prior = c(1,1,2),
                                      misClassCost = NULL,
                                      response = factor(LETTERS[c(1,1,1,2,2,3)]),
                                      insideNode = FALSE)
  expect_equal(length(unique(test2$prior / c(2,3,12))), 1)

  # Subset the classes
  priorAndMisClassCost <- list(prior = structure(c(A = 1, B = 2, C = 3), class = "table", dim = 3L, dimnames = list(
    c("A", "B", "C"))), misClassCost = structure(c(0, 1, 2, 3, 0, 4, 5, 6, 0), dim = c(3L, 3L),
                                                 dimnames = list(c("A", "B", "C"), c("A", "B", "C"))))
  test3 <- updatePriorAndMisClassCost(prior = priorAndMisClassCost$prior,
                                      misClassCost = priorAndMisClassCost$misClassCost,
                                      response = factor(LETTERS[c(1,1,1,3)]),
                                      insideNode = TRUE)
  expect3 <- list(prior = structure(c(A = 0.5, C = 0.5), class = "table", dim = 2L, dimnames = list(
    c("A", "C"))), misClassCost = structure(c(0, 2, 5, 0), dim = c(2L,
                                                                   2L), dimnames = list(c("A", "C"), c("A", "C"))))
  expect_equal(test3, expect3)
})
