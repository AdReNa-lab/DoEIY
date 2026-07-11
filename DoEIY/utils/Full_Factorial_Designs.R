Full_Factorial_Designs <- function(Factor_levels) {
  # Function to create a full factorial design matrix
  # Input: Factor_levels - a numeric vector where each element specifies
  #                        the number of levels for a factor
  # Example: Factor_levels = c(2, 3, 2, 4) means:
  #   - Factor 1 has 2 levels
  #   - Factor 2 has 3 levels
  #   - Factor 3 has 2 levels
  #   - Factor 4 has 4 levels
  
  num_factors = length(Factor_levels)
  
  # Initialize an empty list to store levels for each factor
  All_Levels = list()
  
  # Loop through each factor to create a sequence of levels
  for (i in 1:num_factors) {
    num_levels = Factor_levels[i]
    All_Levels[[length(All_Levels)+1]] = seq(1, num_levels,1)
  }
  
  # Generate the full factorial design by taking all combinations of levels
  Design = expand.grid(All_Levels)
 
  return(Design) 
}
