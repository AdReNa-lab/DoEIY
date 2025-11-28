D_Optimal_Augment <- function(existing_design,
                              levels,
                              components,
                              factor_types,
                              nruns,
                              randomize) {
  
  #
  # 1. Normalise factor names
  #

  var_names <- trimws(unique(all.vars(as.formula(paste("~", paste(components, collapse = " + "))))))
  
  # Keep only relevant factor columns
  colnames(existing_design) <- trimws(colnames(existing_design))
  existing_design <- existing_design[, intersect(var_names, colnames(existing_design)), drop = FALSE]
  
  #
  # 2. Identify polynomial degree of each factor
  #

  factor_degree <- function(fname) {
    deg <- 1
    for (term in components) {
      parts <- strsplit(term, ":")[[1]]
      count <- sum(parts == fname)
      if (count > deg) deg <- count
    }
    deg
  }
  
  #
  # 3. Build full candidate set for safe model parameter counting
  #

  candidate_list <- lapply(var_names, function(f) {

    ftype <- factor_types[[f]]
    levs  <- levels[[f]]
    deg   <- factor_degree(f)
    
    if (ftype == "Continuous") {
      seq(-1, 1, length.out = deg + 1)
      
    } else if (ftype == "Discrete") {
      k <- length(levs)
      sort(as.numeric(levs))
      
    } else if (ftype == "Categorical") {
      factor(levs, levels = levs)
      
    } else {
      stop("Unknown factor type for ", f)
    }
  })
  
  candidate_set <- do.call(expand.grid, candidate_list)
  colnames(candidate_set) <- var_names
  #
  # 4. Fix component strings to ensure quadratic terms enter the formula
  #

  fix_component <- function(term) {
    parts <- strsplit(term, ":")[[1]]
    uniq  <- unique(parts)
    
    if (length(uniq) == 1 && length(parts) > 1) {
      fname <- uniq[1]
      deg   <- length(parts)
      return(paste0("I(", fname, "^", deg, ")"))
    }
    
    term
  }
  
  components_fixed <- sapply(components, fix_component)
  formula <- as.formula(paste("~", paste(components_fixed, collapse = " + ")))
  
  #
  # 5. Compute full model matrix to get number of parameters p
  #

  X_full <- model.matrix(formula, data = candidate_set)
  p <- ncol(X_full)
  
  #
  # 6. Pad the existing design so that model.matrix never fails
  #

  typed_existing <- data.frame(
    row_id = seq_len(nrow(existing_design)),
    stringsAsFactors = FALSE
  )
  
  for (f in var_names) {
    ftype <- factor_types[[f]]
    levs  <- levels[[f]]
    
    if (f %in% colnames(existing_design)) {
      
      col <- existing_design[[f]]
      
      if (ftype %in% c("Continuous", "Discrete")) {
        typed_existing[[f]] <- as.numeric(col)
        
      } else if (ftype == "Categorical") {
        typed_existing[[f]] <- factor(as.character(col), levels = levs)
        
      }
      
    } else {
      
      if (ftype %in% c("Continuous", "Discrete")) {
        typed_existing[[f]] <- rep(0, nrow(existing_design))
      } else if (ftype == "Categorical") {
        typed_existing[[f]] <- factor(rep(levs[1], nrow(existing_design)), levels = levs)
      }
    }
  }

  scale_existing <- function(df, factor_types, levels) {
    
    out <- df
    
    # intersect() ensures we only process factors that actually exist in the data frame
    common_factors <- intersect(names(out), names(factor_types))
    
    for (f in common_factors) {
      ftype <- factor_types[[f]]
      levs  <- levels[[f]]
      col   <- out[[f]]
      
      if (ftype %in% c("Continuous", "Discrete")) {
        
        x <- suppressWarnings(as.numeric(col))
        xmin <- min(x, na.rm = TRUE)
        xmax <- max(x, na.rm = TRUE)
        
        if (xmin == xmax) {
          out[[f]] <- rep(0, length(x))
        } else {
          out[[f]] <- 2 * (x - xmin) / (xmax - xmin) - 1
        }
        
      } else if (ftype == "Categorical") {
        
        # Only map levels that exist in user-supplied levs
        out[[f]] <- factor(
          match(as.character(col), levs),
          levels = seq_len(length(levs))
        )
      }
    }
    
    out
  }
  
  typed_existing <- scale_existing(typed_existing, factor_types, levels)
  
  typed_existing$row_id <- NULL
  
  #
  # 7. Rank of existing design
  #

  X_exist <- model.matrix(formula, data = typed_existing)
  rank_existing <- qr(X_exist)$rank
  
  #
  # 8. Required number of additional runs
  #

  required_new_runs <- max(0, p - rank_existing)
  
  if (nruns < required_new_runs) {
    stop("Model requires ", required_new_runs, " additional runs.")
  }
  
  #
  # 9. Prepare augmented candidate pool
  #

  n_existing <- nrow(typed_existing)
  augmented_candidates <- rbind(typed_existing, candidate_set)
  existing_rows <- seq_len(n_existing)
  
  #
  # 10. Run D optimal augmentation
  #

  max_attempts <- 100
  success <- FALSE
  
  all_n_runs = nruns + n_existing
  
  for (attempt in seq_len(max_attempts)) {
    # Different RNG each time
    
    #print(sample.int(.Machine$integer.max, 1))
    try({
      Design <- optFederov(
        formula,
        data    = augmented_candidates,
        nTrials = all_n_runs,
        augment = TRUE,
        rows    = existing_rows,
        args = TRUE
      )$design

      success <- TRUE
      break
    }, silent = TRUE)
  }

  if (!success) {
    showNotification(
      "The optimisation failed to find a nonsingular D-optimal design after multiple attempts.",
      type = "error",
      duration = 10
    )
    return(NULL)   
  }

  #
  # 11. Reorder existing runs first, then new runs
  #

  is_existing <- apply(Design, 1, function(rD) {
    any(apply(typed_existing, 1, function(rE) all(rD == rE)))
  })

  existing_found <- Design[is_existing, , drop = FALSE]

  order_existing <- sapply(
    seq_len(nrow(typed_existing)),
    function(i) which(apply(existing_found, 1, function(r) all(r == typed_existing[i, ])))
  )

  order_existing <- unlist(order_existing, use.names = FALSE)
  existing_reordered <- existing_found[order_existing, , drop = FALSE]
  
  new_rows <- Design[!is_existing, , drop = FALSE]

  if (randomize && nrow(new_rows) > 1) {
    new_rows <- new_rows[sample(nrow(new_rows)), , drop = FALSE]
  }

  final <- rbind(existing_reordered, new_rows)
  rownames(final) <- NULL
  
  unscale_from_levels <- function(df, factor_types, levels) {
    out <- df
    
    for (f in names(out)) {
      ftype <- factor_types[[f]]
      levs  <- levels[[f]]
      
      if (ftype == "Continuous") {
        
        # Convert coded [-1,1] back to original range
        x <- as.numeric(out[[f]])
        lo <- min(as.numeric(levs))
        hi <- max(as.numeric(levs))
        
        if (lo == hi) {
          out[[f]] <- rep(lo, length(x))
        } else {
          out[[f]] <- ((x + 1) / 2) * (hi - lo) + lo
        }
        
      } else if (ftype == "Discrete") {
        
        # Step 1: unscale like continuous
        x <- as.numeric(out[[f]])
        orig_vals <- as.numeric(levs)
        lo <- min(orig_vals)
        hi <- max(orig_vals)
        
        if (lo == hi) {
          mapped <- rep(lo, length(x))
        } else {
          cont <- ((x + 1) / 2) * (hi - lo) + lo
          
          # Step 2: snap to nearest discrete level
          mapped <- sapply(cont, function(v) {
            orig_vals[which.min(abs(orig_vals - v))]
          })
        }
        out[[f]] <- mapped
        
      } else if (ftype == "Categorical") {
        
        # coded levels are 1..k
        coded <- as.numeric(out[[f]])
        out[[f]] <- factor(levs[coded], levels = levs)
      }
    }
    
    out
  }
  
  final_unscaled <- unscale_from_levels(final, factor_types, levels)
  
  added_factors <- setdiff(var_names, colnames(existing_design))
  
  # Blank added factors for original rows
  if (length(added_factors) > 0) {
    for (f in added_factors) {
      final_unscaled[seq_len(n_existing), f] <- NA
    }
  }
  
  return(final_unscaled)
}
