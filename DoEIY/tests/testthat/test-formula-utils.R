test_that("classify_terms works correctly", {
  expect_equal(classify_terms("A"), 1)
  expect_equal(classify_terms("A:B"), 2)
  expect_equal(classify_terms("A:B:C"), 3)
})

test_that("factor_degree finds correct maximum factor occurrence", {
  components <- c("x1", "x2", "x1:x1", "x1:x2", "x1:x1:x1")
  expect_equal(factor_degree("x1", components), 3)
  expect_equal(factor_degree("x2", components), 1)
  expect_equal(factor_degree("x3", components), 1) # Not present, default is 1
})

test_that("fix_component formats formulas correctly", {
  expect_equal(fix_component("x1"), "x1")
  expect_equal(fix_component("x1:x1"), "I(x1^2)")
  expect_equal(fix_component("x1:x1:x1"), "I(x1^3)")
  expect_equal(fix_component("x1:x2"), "x1:x2")
})

test_that("term_params computes parameters count based on factor type", {
  factor_info <- list(
    x1 = list(name = "x1", type = "Continuous", levels = c("-1", "1")),
    x2 = list(name = "x2", type = "Categorical", levels = c("A", "B", "C")),
    x3 = list(name = "x3", type = "Discrete", levels = c("1", "2", "3"))
  )
  
  # Continuous main effect requires 1 param
  expect_equal(term_params("x1", factor_info), 1)
  
  # Categorical (3 levels) requires k-1 = 2 params
  expect_equal(term_params("x2", factor_info), 2)
  
  # Discrete (3 levels) is handled similarly to Categorical in term_params (levels-1) = 2 params
  expect_equal(term_params("x3", factor_info), 2)
  
  # Interaction A:B requires 1 * 2 = 2 params
  expect_equal(term_params("x1:x2", factor_info), 2)
})

test_that("validate_num_runs validates input runs properly", {
  # Valid integer within bounds
  res <- validate_num_runs(10, 5, 15)
  expect_true(res$valid)
  expect_null(res$message)
  
  # Normalize character input
  res <- validate_num_runs("10", 5, 15)
  expect_true(res$valid)
  
  # Out of bounds (too low)
  res <- validate_num_runs(4, 5, 15)
  expect_false(res$valid)
  expect_match(res$message, "must be within the minumum and maximum")
  
  # Out of bounds (too high)
  res <- validate_num_runs(16, 5, 15)
  expect_false(res$valid)
  
  # Decimals
  res <- validate_num_runs(10.5, 5, 15)
  expect_false(res$valid)
  expect_match(res$message, "must be an integer")
  
  # Zero or negative
  res <- validate_num_runs(0, 0, 15)
  expect_false(res$valid)
  expect_match(res$message, "must be greater than zero")
  
  # Non-numeric character
  res <- validate_num_runs("invalid", 5, 15)
  expect_false(res$valid)
  expect_match(res$message, "must be a single numeric integer")
})
