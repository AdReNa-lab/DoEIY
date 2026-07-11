test_that("Full_Factorial_Designs works with simple inputs", {
  # 2 factors with 2 and 3 levels
  levels <- c(2, 3)
  design <- Full_Factorial_Designs(levels)
  
  expect_s3_class(design, "data.frame")
  expect_equal(ncol(design), 2)
  expect_equal(nrow(design), 6) # 2 * 3 = 6
  
  # Verify specific values are present
  expect_equal(design[[1]], c(1, 2, 1, 2, 1, 2))
  expect_equal(design[[2]], c(1, 1, 2, 2, 3, 3))
})

test_that("Full_Factorial_Designs works with a single factor", {
  design <- Full_Factorial_Designs(c(5))
  expect_equal(nrow(design), 5)
  expect_equal(ncol(design), 1)
  expect_equal(design[[1]], 1:5)
})
