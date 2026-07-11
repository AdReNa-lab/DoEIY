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

transform_term <- function(term) {
  # Transforms term formatting: quadratics (from A * A to I(A^2)) and interactions (from A * B to A:B)
  if (grepl("\\*", term)) {
    vars <- unlist(strsplit(term, "\\s*\\*\\s*"))
    if (length(vars) == 2 && vars[1] == vars[2]) {
      # Quadratic term like "A * A"
      paste0("I(", vars[1], "^2)")
    } else {
      # Interaction term like "A * B"
      paste0(vars[1], ":", vars[2])
    }
  } else {
    term
  }
}

update_factor_names <- function(new_columns, design_matrix) {
  # Cleans up column names after model analysis
  original_columns <- colnames(design_matrix)
  # Initialize a named list to track renamed columns
  renamed_columns <- list()

  # Step 1: Rename columns if they match the pattern and have 2 unique values
  for (col in new_columns) {
    for (original_col in original_columns) {
      # Create a pattern to match the original column name followed by a number
      pattern <- paste0("^", original_col, "[0-9]+$")

      # Check if the current column name matches this pattern
      if (grepl(pattern, col)) {
        # Find out how many unique values are in the corresponding column in design_matrix
        unique_values <- unique(design_matrix[[original_col]])
        num_unique_values <- length(unique_values)

        # If there are exactly 2 unique values, rename the column in new_columns
        if (num_unique_values == 2) {
          new_name <- original_col
          new_columns[new_columns == col] <- new_name
          renamed_columns[[col]] <- new_name
        }
      }
    }
  }

  # Step 2: Update interaction terms in new_columns based on renamed columns
  for (col in new_columns) {
    for (old_name in names(renamed_columns)) {
      # Replace interaction terms in the format "old_name:other" with "new_name:other"
      # but ensure only whole-word matches are replaced
      new_name <- renamed_columns[[old_name]]
      updated_col <- gsub(paste0("\\b", old_name, "\\b"), new_name, col)

      # Update the column name if it was changed
      if (updated_col != col) {
        new_columns[new_columns == col] <- updated_col
      }
    }
  }

  # Step 3: Loop through each column name in new_columns and update interaction notation
  for (col in new_columns) {
    if (grepl(":", col)) {
      # Split the column name by ":"
      terms <- unlist(strsplit(col, ":"))

      # Check if all possible combinations of terms match any name in original_columns
      valid_combination <- TRUE
      for (i in seq_along(terms)) {
        combined_term <- terms[i]

        # Start combining terms one by one to check for multi-part names
        for (j in (i + 1):length(terms)) {
          combined_term <- paste0(combined_term, ":", terms[j])

          # Check if this combined term exists in original_columns
          if (combined_term %in% original_columns) {
            terms <- c(terms[1:(i - 1)], combined_term, terms[(j + 1):length(terms)])
            break
          }
        }

        # If any part of the name isn't valid, stop the process
        if (!(terms[i] %in% original_columns)) {
          valid_combination <- FALSE
          break
        }
      }

      # If the combination is valid, replace ":" with "*" in the column name
      if (valid_combination) {
        new_name <- paste(terms, collapse = "*")
        new_columns[new_columns == col] <- new_name
      }
    }
  }

  for (col in new_columns) {
    if (grepl("\\.", col)) {
      # Split the column name by "."
      terms <- unlist(strsplit(col, "\\."))

      # Check if all possible combinations of terms match any name in original_columns
      valid_combination <- TRUE
      for (i in seq_along(terms)) {
        combined_term <- terms[i]

        # Start combining terms one by one to check for multi-part names
        for (j in (i + 1):length(terms)) {
          combined_term <- paste0(combined_term, ".", terms[j])

          # Check if this combined term exists in original_columns
          if (combined_term %in% original_columns) {
            terms <- c(terms[1:(i - 1)], combined_term, terms[(j + 1):length(terms)])
            break
          }
        }

        # If any part of the name isn't valid, stop the process
        if (!(terms[i] %in% original_columns)) {
          valid_combination <- FALSE
          break
        }
      }

      # If the combination is valid, replace "." with "*" in the column name
      if (valid_combination) {
        new_name <- paste(terms, collapse = "*")
        new_columns[new_columns == col] <- new_name
      }
    }
  }
  # Return the updated column names
  new_columns
}

