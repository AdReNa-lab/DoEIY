Latin_Hypercube_Designs <- function(num_factors, num_runs) {
  
  # Generate the LHS design using the improvedLHS method for better space-filling
  Design <- improvedLHS(n = num_runs, k = num_factors)
  
  # Convert the design to a data frame and label the columns
  Design <- as.data.frame(Design)
  colnames(Design) <- paste0("Var", seq_len(num_factors))
  
  return(Design)
}