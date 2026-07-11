Fractional_Factorial_Designs <- function(num_runs, num_factors, num_blocks) {
  # All inputs are integers
  
  # Build the fractional factorial design with FrF2
  Design <- as.data.frame(FrF2(nruns = num_runs, nfactors = num_factors, blocks = num_blocks))
  
  # Convert any factor columns to numeric; keep non-factors as is
  Design[] <- lapply(Design , function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  
  
  # Drop Blocks column if present (renaming handled upstream as you noted)
  if ("Blocks" %in% colnames(Design)) {
    design_matrix <- subset(Design, select = -c(Blocks))
  }  else {design_matrix <- Design}
  
  
  p  <- ncol(design_matrix)
  n  <- nrow(design_matrix)
  cn <- colnames(design_matrix)
  
  
  # Build 2-factor interactions (for p < 2)
  if (p >= 2L) {
    two_idx <- utils::combn(p, 2, simplify = TRUE)
    two_factor_interactions <- sapply(seq_len(ncol(two_idx)), function(k) {
      i <- two_idx[1L, k]; j <- two_idx[2L, k]
      design_matrix[[i]] * design_matrix[[j]]
    })
    if (is.null(dim(two_factor_interactions))) two_factor_interactions <- matrix(two_factor_interactions, ncol = 1L)
    colnames(two_factor_interactions) <- apply(two_idx, 2, function(v) paste(cn[v], collapse = ":"))
  } else {
    two_factor_interactions <- matrix(numeric(0), nrow = n, ncol = 0)
  }
  
  # Build 3-factor interactions (for p < 3)
  if (p >= 3L) {
    three_idx <- utils::combn(p, 3, simplify = TRUE)
    three_factor_interactions <- sapply(seq_len(ncol(three_idx)), function(k) {
      i <- three_idx[1L, k]; j <- three_idx[2L, k]; l <- three_idx[3L, k]
      design_matrix[[i]] * design_matrix[[j]] * design_matrix[[l]]
    })
    if (is.null(dim(three_factor_interactions))) three_factor_interactions <- matrix(three_factor_interactions, ncol = 1L)
    colnames(three_factor_interactions) <- apply(three_idx, 2, function(v) paste(cn[v], collapse = ":"))
  } else {
    three_factor_interactions <- matrix(numeric(0), nrow = n, ncol = 0)
  }
  

  # Helper Function - check if aliased
  is_aliased <- function(a, b, tol = 1e-12) {
    # keep only finite pairs
    ok <- is.finite(a) & is.finite(b)
    a <- a[ok]
    b <- b[ok]
    if (!length(a)) return(FALSE)
    
    # Get SD to see if columns are essentially constant (accounts for floating-point noise)
    sa <- stats::sd(a)
    sb <- stats::sd(b)
    
    # if both (near) constant they are consider aliased 
    if (sa <= tol && sb <= tol) return(TRUE)
    
    # if one constant and the other not they're not aliased under 
    if (sa <= tol || sb <= tol) return(FALSE)
    
    # Check for perfect linear relationship (with or without intercept) 
    # If two columns are equal or multiples of each other, r should be 1
    r <- suppressWarnings(stats::cor(a, b, use = "complete.obs"))
    # accounts for floating-point noise
    is.finite(r) && (abs(r) >= 1 - tol)
  }
  
  
  # Resolution checks
  # Main factor aliasing
  main_aliases <- 0L
  if (p >= 2L) {
    for (i in 1:(p - 1L)) {
      for (j in (i + 1L):p) {
        if (is_aliased(design_matrix[, i], design_matrix[, j])) { main_aliases <- 1L; break }
      }
      if (main_aliases) break
    }
  }
  
  # Check if main factors are aliased with 2-factor interactions
  main_secondary_aliases <- 0L
  if (!main_aliases && ncol(two_factor_interactions) > 0L) {
    for (i in 1:p) {
      for (j in 1:ncol(two_factor_interactions)) {
        if (is_aliased(design_matrix[, i], two_factor_interactions[, j])) { main_secondary_aliases <- 1L; break }
      }
      if (main_secondary_aliases) break
    }
  }
  
  # Check if 2-factor interactions are aliased with each other
  secondary_aliases <- 0L
  if (!main_aliases && !main_secondary_aliases && ncol(two_factor_interactions) >= 2L) {
    for (i in 1:(ncol(two_factor_interactions) - 1L)) {
      for (j in (i + 1L):ncol(two_factor_interactions)) {
        if (is_aliased(two_factor_interactions[, i], two_factor_interactions[, j])) { secondary_aliases <- 1L; break }
      }
      if (secondary_aliases) break
    }
  }
  
  # Check if 2-factor interactions are aliased with 3-factor interactions
  secondary_tertiary_aliases <- 0L
  if (!main_aliases && !main_secondary_aliases && !secondary_aliases &&
      ncol(two_factor_interactions) > 0L && ncol(three_factor_interactions) > 0L) {
    for (i in 1:ncol(two_factor_interactions)) {
      for (j in 1:ncol(three_factor_interactions)) {
        if (is_aliased(two_factor_interactions[, i], three_factor_interactions[, j])) { secondary_tertiary_aliases <- 1L; break }
      }
      if (secondary_tertiary_aliases) break
    }
  }
  

  # Determine the resolution of the design based on the aliasing
  if (main_aliases == 1) {
    resolution = "Resolution II"
  } else if (main_secondary_aliases == 1) {
    resolution = "Resolution III"
  } else if (secondary_aliases == 1) {
    resolution = "Resolution IV"
  } else if (secondary_tertiary_aliases == 1) {
    resolution = "Resolution V"
  } else {
    resolution = "Resolution VI or higher"
  }

  
  return(list(Design = Design, Resolution = resolution))
}