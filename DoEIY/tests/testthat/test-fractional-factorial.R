test_that("Fractional_Factorial_Designs generates valid fractional factorial designs", {
  library(FrF2)
  
  # 8 runs, 4 factors, no blocking
  res <- Fractional_Factorial_Designs(8, 4, 1)
  
  expect_type(res, "list")
  expect_named(res, c("Design", "Resolution"))
  
  design <- res$Design
  expect_s3_class(design, "data.frame")
  expect_equal(ncol(design), 4)
  expect_equal(nrow(design), 8)
  
  # Verify resolution type
  expect_type(res$Resolution, "character")
  expect_true(grepl("Resolution", res$Resolution))
})
