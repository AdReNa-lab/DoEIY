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

test_that("transform_term formats polynomial and interaction terms correctly", {
  expect_equal(transform_term("A * A"), "I(A^2)")
  expect_equal(transform_term("A * B"), "A:B")
  expect_equal(transform_term("A"), "A")
})

test_that("update_factor_names cleans and formats ANOVA column names correctly", {
  design_matrix <- data.frame(
    A = c(-1, 1, -1, 1),
    B = factor(c("Low", "Medium", "High", "Low")),
    C = c(10, 20, 30, 40)
  )
  
  # A has exactly 2 unique values, so A1 -> A
  # B has 3 unique values, so B1 is NOT renamed to B
  # A1:B should update to A:B, then since A and B are original columns, A:B -> A*B
  # A1.B should update to A.B, then since A and B are original columns, A.B -> A*B
  
  new_cols <- c("A1", "B1", "A1:B", "A1.B", "C")
  updated_cols <- update_factor_names(new_cols, design_matrix)
  
  expect_equal(updated_cols[1], "A")
  expect_equal(updated_cols[2], "B1")
  expect_equal(updated_cols[3], "A*B")
  expect_equal(updated_cols[4], "A*B")
  expect_equal(updated_cols[5], "C")
})

