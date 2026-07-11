Box_Behnken_Designs <- function(n) {
  # Design Matrix Skeleton
  N03_Design_Matrix <- matrix(c(
    1, 1, 0,
    1, 0, 1,
    0, 1, 1
  ), ncol = 3, byrow = TRUE)
  N03_Blocks <- c(3)          # one block containing all 3 rows
  N03_CPs    <- c(3)          # 3 center points total
  
  
  N04_Design_Matrix <- matrix(c(
    1, 1, 0, 0,
    0, 0, 1, 1,
    1, 0, 0, 1,
    0, 1, 1, 0,
    1, 0, 1, 0,
    0, 1, 0, 1
  ), ncol = 4, byrow = TRUE)
  N04_Blocks <- c(2, 2, 2)    # 3 blocks of 2 rows each
  N04_CPs    <- c(1, 1, 1)    # 1 center point per block
  
  
  N05_Design_Matrix <- matrix(c(
    1, 1, 0, 0, 0,
    0, 0, 1, 1, 0,
    0, 1, 0, 0, 1,
    1, 0, 1, 0, 0,
    0, 0, 0, 1, 1,
    0, 1, 1, 0, 0,
    1, 0, 0, 1, 0,
    0, 0, 1, 0, 1,
    1, 0, 0, 0, 1, 
    0, 1, 0, 1, 0
  ), ncol = 5, byrow = TRUE)
  N05_Blocks <- c(5, 5)       # 2 blocks of 5 rows each
  N05_CPs    <- c(3, 3)       # 3 center points per block
  
  
  N06_Design_Matrix <- matrix(c(
    1, 1, 0, 1, 0, 0,
    0, 1, 1, 0, 1, 0,
    0, 0, 1, 1, 0, 1,
    1, 0, 0, 1, 1, 0,
    0, 1, 0, 0, 1, 1,
    1, 0, 1, 0, 0, 1
  ), ncol = 6, byrow = TRUE)
  N06_Blocks <- c(6)          # FIXED: one block containing all 6 rows
  N06_CPs    <- c(6)          # 6 center points total
  
  
  N07_Design_Matrix <- matrix(c(
    0, 0, 0, 1, 1, 1, 0,
    1, 0, 0, 0, 0, 1, 1,
    0, 1, 0, 0, 1, 0, 1,
    1, 1, 0, 1, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 1,
    1, 0, 1, 0, 1, 0, 0,
    0, 1, 1, 0, 0, 1, 0
  ), ncol = 7, byrow = TRUE)
  N07_Blocks <- c(7)          # FIXED: one block containing all 7 rows
  N07_CPs    <- c(6)          # 6 center points total
  

  
  pick <- switch(
    as.character(n),
    "3" = list(M = N03_Design_Matrix, Blocks = N03_Blocks, CPs = N03_CPs),
    "4" = list(M = N04_Design_Matrix, Blocks = N04_Blocks, CPs = N04_CPs),
    "5" = list(M = N05_Design_Matrix, Blocks = N05_Blocks, CPs = N05_CPs),
    "6" = list(M = N06_Design_Matrix, Blocks = N06_Blocks, CPs = N06_CPs),
    "7" = list(M = N07_Design_Matrix, Blocks = N07_Blocks, CPs = N07_CPs),
    NULL
  )
  if (is.null(pick)) {
    stop("")
  }
  
  Design_Matrix <- pick$M
  Blocks <- pick$Blocks
  CPs    <- pick$CPs
  
  
  if (sum(Blocks) != nrow(Design_Matrix)) {
    stop("sum(Blocks) must equal nrow(Design_Matrix).")
  }
  if (length(CPs) != length(Blocks)) {
    stop("length(CPs) must equal length(Blocks).")
  }
  
  colnames(Design_Matrix) <- paste0("X", seq_len(ncol(Design_Matrix)))
  row_ptr <- 0L
  out <- NULL
  
  for (i in seq_along(Blocks)) {
    rows_in_block <- Blocks[i]
    block_runs <- NULL
    
    for (j in seq_len(rows_in_block)) {
      row_idx <- row_ptr + j
      row <- Design_Matrix[row_idx, ]
      
      # Expand active factors -> {-1, +1}, inactive -> {0}
      levels_list <- lapply(as.list(row), function(v) if (v == 1) c(-1, 1) else 0)
      expanded <- expand.grid(levels_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      
      if (is.null(block_runs)) block_runs <- expanded else block_runs <- rbind(block_runs, expanded)
    }
    
    # Append center points for this block
    if (CPs[i] > 0) {
      cps <- as.data.frame(matrix(0, nrow = CPs[i], ncol = ncol(Design_Matrix)))
      names(cps) <- colnames(Design_Matrix)
      block_runs <- rbind(block_runs, cps)
    }
    
    block_runs$Block <- i
    if (is.null(out)) out <- block_runs else out <- rbind(out, block_runs)
    
    row_ptr <- row_ptr + rows_in_block
  }
  
  rownames(out) <- NULL
  out
}