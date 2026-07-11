# Test runner for DoEIY design functions
library(testthat)

# Source all the design generators
source("Box_Behnken_Designs.R")
source("Full_Factorial_Designs.R")
source("Plackett_Burman_Designs.R")
source("Fractional_Factorial_Designs.R")
source("Central_Composite_Designs.R")
source("Latin_Hypercube_Designs.R")
source("D_Optimal_Design.R")
source("D_Optimal_Augment.R")

# Run all tests in the tests/testthat directory
test_dir("tests/testthat", reporter = "summary")
