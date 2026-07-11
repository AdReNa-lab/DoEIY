test_that("Plackett_Burman_Designs creates valid designs for supported factors", {
  # Plackett-Burman designs are constructed with specific run sizes (multiples of 4)
  # For factors 4 to 7, it constructs an 8-run design (+ 1 row of all -1s = 8 rows total)
  design_7 <- Plackett_Burman_Designs(7)
  expect_s3_class(design_7, "data.frame")
  expect_equal(ncol(design_7), 7)
  expect_equal(nrow(design_7), 8)
  
  # For 8 factors, it should use the 12-run design (+ 1 row of all -1s = 12 rows)
  design_8 <- Plackett_Burman_Designs(8)
  expect_equal(ncol(design_8), 8)
  expect_equal(nrow(design_8), 12)
  
  # For 12 factors, it uses 16 runs
  design_12 <- Plackett_Burman_Designs(12)
  expect_equal(ncol(design_12), 12)
  expect_equal(nrow(design_12), 16)
})

test_that("Plackett_Burman_Designs handles limits and invalid factors", {
  expect_error(Plackett_Burman_Designs(3))
  expect_error(Plackett_Burman_Designs(24))
})
