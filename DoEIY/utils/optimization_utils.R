# Optimization and evaluation utilities for DoEIY

objective_value <- function(pred, opt_type, target = NULL) {
  # Standardizes predictions to a minimization objective
  if (opt_type == "Maximize") {
    return(-pred)
  }
  if (opt_type == "Minimize") {
    return(pred)
  }
  if (opt_type == "Find Target") {
    if (is.null(target) || is.na(target)) {
      stop("Target response is required.")
    }
    return((pred - target)^2)
  }
  return(pred)
}
