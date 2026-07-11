test_that("Central_Composite_Designs creates valid designs", {
  # We need the rsm package loaded for these functions since they call rsm's ccd
  library(rsm)
  
  # Circumscribed
  design_cc <- Central_Composite_Designs(3, "Circumscribed")
  expect_s3_class(design_cc, "data.frame")
  expect_equal(ncol(design_cc), 3)
  
  # Inscribed
  design_in <- Central_Composite_Designs(3, "Inscribed")
  expect_s3_class(design_in, "data.frame")
  expect_equal(ncol(design_in), 3)
  
  # Face Centered
  design_fc <- Central_Composite_Designs(3, "Face Centered")
  expect_s3_class(design_fc, "data.frame")
  expect_equal(ncol(design_fc), 3)
  
  # Confirm face-centered levels are strictly within [-1, 1]
  expect_true(all(design_fc >= -1 & design_fc <= 1))
})

test_that("Central_Composite_Designs throws error on unknown type", {
  expect_error(Central_Composite_Designs(3, "UnknownType"))
})
