test_that("Box_Behnken_Designs generates valid designs for supported factors", {
  # Test with 3 factors
  design_3 <- Box_Behnken_Designs(3)
  expect_s3_class(design_3, "data.frame")
  expect_true("Block" %in% colnames(design_3))
  # Standard 3-factor Box-Behnken design usually has 12 factorial points + 3 center points = 15 runs
  expect_equal(nrow(design_3), 15)
  expect_equal(ncol(design_3), 4) # X1, X2, X3, Block
  
  # Test with 4 factors (usually 24 factorial points + 3 center points = 27 runs)
  design_4 <- Box_Behnken_Designs(4)
  expect_equal(nrow(design_4), 27)
  
  # Test for factors 5, 6, 7
  expect_no_error(Box_Behnken_Designs(5))
  expect_no_error(Box_Behnken_Designs(6))
  expect_no_error(Box_Behnken_Designs(7))
})

test_that("Box_Behnken_Designs handles invalid factors", {
  expect_error(Box_Behnken_Designs(2))
  expect_error(Box_Behnken_Designs(8))
  expect_error(Box_Behnken_Designs("invalid"))
})
