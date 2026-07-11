Central_Composite_Designs <- function(num_factors, ccd_type) {
  # Generates a Central Composite Design (CCD) matrix
  
  # Define the CCD design based on the type
  Design <- switch(
    ccd_type,
    "Circumscribed" = ccd(num_factors, alpha = "orthogonal", randomize = FALSE),
    "Inscribed"     = ccd(num_factors, alpha = "orthogonal", randomize = FALSE, inscribed = TRUE),
    "Face Centered" = ccd(num_factors, alpha = "face",       randomize = FALSE),
    {
      stop("Unknown ccd_type. Use 'Circumscribed', 'Inscribed', or 'Face Centered'.")
    }
  )
  Design <- as.data.frame(unclass(Design), stringsAsFactors = FALSE)
  
  # Drop the blocking column, not currently supported (also drop other unnecessary columns)
  drop_cols <- intersect(c("run.order", "std.order", "Block"), colnames(Design))
  if (length(drop_cols)) {
    Design <- Design[, !colnames(Design) %in% drop_cols, drop = FALSE]
  }
  
  colnames(Design) <- paste0("Var", seq_len(ncol(Design)))
  
  return(Design)
}