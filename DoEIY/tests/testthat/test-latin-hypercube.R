test_that("Latin_Hypercube_Designs generates valid designs", {
  library(lhs)
  
  num_factors <- 3
  num_runs <- 10
  
  design <- Latin_Hypercube_Designs(num_factors, num_runs)
  
  expect_s3_class(design, "data.frame")
  expect_equal(ncol(design), num_factors)
  expect_equal(nrow(design), num_runs)
  
  # Confirm space filling values are strictly within [0, 1]
  expect_true(all(design >= 0 & design <= 1))
  
  # Column names should be Var1, Var2, Var3
  expect_equal(colnames(design), c("Var1", "Var2", "Var3"))
})
