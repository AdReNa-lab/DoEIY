test_that("get_design_description returns correct HTML for valid designs", {
  # Valid designs
  desc_pb <- get_design_description("Plackett-Burman")
  expect_true(grepl("Plackett-Burman Designs", desc_pb))
  
  desc_ff <- get_design_description("Full Factorial")
  expect_true(grepl("Full Factorial Designs", desc_ff))
  
  desc_do <- get_design_description("D-Optimal")
  expect_true(grepl("D-Optimal Designs", desc_do))
})

test_that("get_design_description returns empty string for invalid design", {
  desc_invalid <- get_design_description("Non-Existent Design")
  expect_equal(desc_invalid, "")
  
  desc_empty <- get_design_description("")
  expect_equal(desc_empty, "")
})
