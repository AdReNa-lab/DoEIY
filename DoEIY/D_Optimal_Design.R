D_Optimal_Designs <- function(levels, components, factor_types, nruns) {

  # Build model formula from supplied terms (components)
  formula <- as.formula(paste("~", paste(components, collapse = " + ")))
  
  # Extract factor names robustly from the formula
  var_names <- unique(all.vars(formula))
  
  num_factors <- length(levels)
  
  factor_degree <- function(fname, components) {
    deg <- 1
    for (term in components) {
      parts <- unlist(strsplit(term, ":"))
      count <- sum(parts == fname)
      if (count > deg) deg <- count
    }
    deg
  }
  
  # Build candidate set
  candidate_list <- lapply(seq_len(num_factors), function(i) {
    fname <- var_names[i]
    ft    <- factor_types[i]
    k     <- levels[i]
    
    deg   <- factor_degree(fname, components)
    
    if (ft == "Continuous") {
      
      if (deg == 1) {
        seq(-1, 1, length.out = 2)
      } else if (deg == 2) {
        seq(-1, 1, length.out = 3)
      } else {
        seq(-1, 1, length.out = deg + 1)
      }
      
    } else if (ft == "Discrete") {
      if (k == 2) c(-1, 1) else seq(-1, 1, length.out = k)
    } else if (ft == "Categorical") {
      factor(seq_len(k))
    } else {
      stop("Unknown factor type '", ft, "'")
    }
  })
      
  candidate_set <- do.call(expand.grid, candidate_list)
  
  # Rename columns to match variables referenced in the formula
  colnames(candidate_set) <- var_names
  
  # Ensure the number of rows in the candidate set exceedes the number of runs selected by the user
  if (nruns > nrow(candidate_set)) {                                
    stop("nruns (", nruns, ") exceeds size of candidate set (", nrow(candidate_set), "). Increase levels or reduce nruns.") 
  }
      
  # Generate the D-optimal design
  Design <- optFederov(formula, data = candidate_set, nTrials = nruns)$design
  rownames(Design) <- NULL
  
  return(Design)
}
