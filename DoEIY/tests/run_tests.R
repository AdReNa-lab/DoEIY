# Test runner for DoEIY design functions
library(testthat)

# Source all the design generators
source("utils/formula_utils.R")
source("utils/optimization_utils.R")
source("utils/design_descriptions.R")
source("utils/Box_Behnken_Designs.R")
source("utils/Full_Factorial_Designs.R")
source("utils/Plackett_Burman_Designs.R")
source("utils/Fractional_Factorial_Designs.R")
source("utils/Central_Composite_Designs.R")
source("utils/Latin_Hypercube_Designs.R")
source("utils/D_Optimal_Design.R")
source("utils/D_Optimal_Augment.R")

# Run all tests in the tests/testthat directory
test_dir("tests/testthat", reporter = "summary")
