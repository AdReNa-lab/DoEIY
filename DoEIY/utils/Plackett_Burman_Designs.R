Plackett_Burman_Designs <- function(n) {
  # Build a two-level Plackett-Burman screening design
  # for n factors using a seeded generator and cyclic permutations.
  
  # stop the code if the number of factors is not supported
  if (n < 4) {
    stop("Designs are not recommended for fewer than 4 factors.")
  }
  if (n > 23) {
    stop("This function currently supports up to 23 factors (24-run design).")
  }
  
  # Choose generator based on number of factors
  if (n >= 4 && n < 8) {
    # 8 runs 
    Generator <- c(+1, +1, +1, -1, 
                   +1, -1, -1)
  } else if (n >= 8 && n < 12) {
    # 12 runs
    Generator <- c(+1, +1, -1, +1, 
                   +1, +1, -1, -1, 
                   -1, +1, -1)
  } else if (n >= 12 && n < 16) {
    # 16 runs
    Generator <- c(+1, +1, +1, +1,
                   -1, +1, -1, +1,
                   +1, -1, -1, +1,
                   -1, -1, -1)
  } else if (n >= 16 && n < 20) {
    # 20 runs
    Generator <- c(+1, +1, -1, -1,
                   +1, +1, +1, +1,
                   -1, +1, -1, +1,
                   -1, -1, -1, -1,
                   +1, +1, -1)
  } else if (n >= 20 && n <= 23) {
    # 24 runs
    Generator <- c(+1, +1, +1, +1,
                   +1, -1, +1, -1,
                   +1, +1, -1, -1,
                   +1, +1, -1, -1,
                   +1, -1, +1, -1,
                   -1, -1, -1)
  } else {
    # Should not be reachable, just in case.
    stop("Error selecting generator for number of factors = ", n)
  }
  
  # Construct an empty matrix first (without the final -1 row)
  Design <- matrix(NA_real_, nrow = length(Generator), ncol = n)
  # First column is the generator
  Design[, 1] <- Generator
  
  # Permute the generator to create the additional columns
  permuted <- Generator
  if (n >= 2) {
    for (j in 2:n) {
      permuted <- c(tail(permuted, 1), head(permuted, -1))  # rotate down by 1
      Design[, j] <- permuted
    }
  }
  
  Design <- rbind(Design, rep(-1, n))
  
  Design = as.data.frame(Design)
  
  return(Design)
}
