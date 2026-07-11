test_that("D_Optimal_Designs creates valid design matrix", {
  library(AlgDesign)
  
  levels <- c(3, 3)
  components <- c("x1", "x2", "x1:x1", "x2:x2", "x1:x2")
  factor_types <- c("Continuous", "Continuous")
  nruns <- 6
  
  design <- D_Optimal_Designs(levels, components, factor_types, nruns)
  
  expect_s3_class(design, "data.frame")
  expect_equal(ncol(design), 2)
  expect_equal(nrow(design), 6)
  expect_equal(colnames(design), c("x1", "x2"))
  
  # Values should be within [-1, 1] for continuous/discrete
  expect_true(all(design >= -1 & design <= 1))
})

test_that("D_Optimal_Augment augments existing design", {
  library(AlgDesign)
  
  # Create a starting design
  existing <- data.frame(
    x1 = c(-1, 1, -1),
    x2 = c(-1, -1, 1)
  )
  
  levels <- list(x1 = c(-1, 0, 1), x2 = c(-1, 0, 1))
  components <- c("x1", "x2", "x1:x1", "x2:x2", "x1:x2")
  factor_types <- list(x1 = "Continuous", x2 = "Continuous")
  nruns <- 3 # add 3 runs
  
  augmented <- D_Optimal_Augment(
    existing_design = existing,
    levels = levels,
    components = components,
    factor_types = factor_types,
    nruns = nruns,
    randomize = FALSE
  )
  
  expect_s3_class(augmented, "data.frame")
  expect_equal(ncol(augmented), 2)
  expect_equal(nrow(augmented), 6) # 3 existing + 3 new = 6
  expect_equal(colnames(augmented), c("x1", "x2"))
  
  # Existing runs should remain unchanged in the first 3 rows
  # Wait, let's verify if D_Optimal_Augment returns original scale or normalized scale
  # D_Optimal_Augment unscales the output, so it should match the original existing values
  expect_equal(augmented$x1[1:3], c(-1, 1, -1))
  expect_equal(augmented$x2[1:3], c(-1, -1, 1))
})

