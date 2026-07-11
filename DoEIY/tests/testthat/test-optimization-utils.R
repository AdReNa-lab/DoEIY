test_that("objective_value transforms prediction correctly based on optimization type", {
  # Maximize -> returns negative of prediction
  expect_equal(objective_value(5, "Maximize"), -5)
  expect_equal(objective_value(-3, "Maximize"), 3)
  
  # Minimize -> returns prediction unchanged
  expect_equal(objective_value(5, "Minimize"), 5)
  expect_equal(objective_value(-3, "Minimize"), -3)
  
  # Find Target -> returns squared difference
  expect_equal(objective_value(8, "Find Target", target = 10), 4)
  expect_equal(objective_value(12, "Find Target", target = 10), 4)
  
  # Other types -> return prediction unchanged
  expect_equal(objective_value(5, "UnknownType"), 5)
})

test_that("objective_value throws error when target is missing for Find Target", {
  expect_error(objective_value(5, "Find Target"))
  expect_error(objective_value(5, "Find Target", target = NULL))
  expect_error(objective_value(5, "Find Target", target = NA))
})
