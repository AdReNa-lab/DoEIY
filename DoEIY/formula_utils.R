# Formula and Model parsing utilities for DoEIY

classify_terms <- function(term) {
  # Classifies the degree/interaction level of a term
  length(unlist(strsplit(term, ":")))
}

factor_degree <- function(fname, components) {
  # Determines the maximum polynomial degree of factor fname across components
  deg <- 1
  for (term in components) {
    parts <- unlist(strsplit(term, ":"))
    count <- sum(parts == fname)
    if (count > deg) deg <- count
  }
  deg
}

fix_component <- function(term) {
  # Standardizes pure polynomial terms into R I() syntax for formulas (e.g. X:X -> I(X^2))
  parts <- strsplit(term, ":")[[1]]
  uniq  <- unique(parts)
  
  if (length(uniq) == 1 && length(parts) > 1) {
    fname <- uniq[1]
    deg   <- length(parts)
    return(paste0("I(", fname, "^", deg, ")"))
  }
  
  term
}

term_params <- function(term, factor_info) {
  # Calculates the number of parameters (degrees of freedom) added by a term
  factors <- unlist(strsplit(term, ":"))
  tab <- table(factors)
  
  prod(sapply(names(tab), function(f) {
    info <- factor_info[[f]]
    if (is.null(info)) {
      stop("Factor information not found for factor: ", f)
    }
    if (info$type == "Continuous") {
      1
    } else {
      length(info$levels) - 1
    }
  }))
}

validate_num_runs <- function(num_runs, min_runs, max_runs) {
  # Validates that the runs input is a valid positive integer within the min/max limits.
  # Returns a list: list(valid = TRUE/FALSE, message = character/NULL)
  
  # Normalize character or factor inputs
  if (is.character(num_runs) || is.factor(num_runs)) {
    num_runs_num <- suppressWarnings(as.numeric(as.character(num_runs)))
  } else {
    num_runs_num <- num_runs
  }
  
  # Check it's numeric
  if (!is.numeric(num_runs_num) || length(num_runs_num) != 1 || is.na(num_runs_num)) {
    return(list(valid = FALSE, message = "The 'Number of Experimental Runs' must be a single numeric integer."))
  }

  # Check it's an integer
  if (abs(num_runs_num - round(num_runs_num)) > .Machine$double.eps^0.5) {
    return(list(valid = FALSE, message = "The 'Number of Experimental Runs' must be an integer (no decimal values)."))
  }

  # Check it's positive
  if (num_runs_num <= 0) {
    return(list(valid = FALSE, message = "The 'Number of Experimental Runs' must be greater than zero."))
  }
  
  # Check bounds
  if (!(num_runs_num >= min_runs && num_runs_num <= max_runs)) {
    return(list(valid = FALSE, message = "The 'Number of Experimental Runs' must be within the minumum and maximum number of runs provided."))
  }
  
  return(list(valid = TRUE, message = NULL))
}
