# ==============================================================================
# M340L: Gaussian Elimination Coding Project
# Raquel Mejia-Trujillo
#
# Instructions:
# Code up row reduction for an m x n matrix.
# Test with 11 x 10 matrix with random (uniformly distributed from 0-1) entries.
# Print the result of row reduction to echelon form.
# ==============================================================================

# Define a row reduction function to simple echelon form
SimpleEchelonElimination = function(x){
  # Define matrix dimensions
  n = ncol(x)
  m = nrow(x)

  # Define starting pivot position
  pivot_row = 1
  pivot_col = 1

  # Set FLOPs counter to zero
  FLOPs = 0

  # ============================================================================
  # Swap two rows
  # FLOPs per use of `SwapRows`: 0
  # ============================================================================
  SwapRows = function(x, original_first_row, original_second_row){
    temp_first_row = x[original_first_row, ]
    x[original_first_row, ] = x[original_second_row, ]
    x[original_second_row, ] = temp_first_row
    return(x)
  }

  # ============================================================================
  # Row elimination
  # FLOPs per row use of `EliminateRows`: 2n + 1
  # ============================================================================
  EliminateRows = function(x, pivot_row, pivot_col, eliminated_row){
    # Calculate scaling factor (1 FLOP per row)
    scaling = x[eliminated_row, pivot_col]/x[pivot_row, pivot_col]

    # Replace the newly eliminated row (2 FLOPs per element)
    x[eliminated_row, ] = x[eliminated_row, ]  - (scaling * x[pivot_row, ])
    return(x)
  }

  # ============================================================================
  # Perform elimination to simple echelon form with partial pivoting
  # ============================================================================
  while (pivot_row <= m & pivot_col <= n) {
    # --------------------------------------------------------------------------
    # Partial pivoting
    # --------------------------------------------------------------------------
    # Current pivot entry
    current_pivot = x[pivot_row, pivot_col]
    # Subset the pivot column from current pivot row to the bottom of the matrix `x`
    subset_below_pivot = x[pivot_row:m, pivot_col]
    # Value of the entry in subset that contains the largest absolute value
    max_entry_pivot_col = max(abs(subset_below_pivot))
    # Index of the row in the original matrix that contains the entry defined above
    new_pivot_row = which(abs(x[, pivot_col]) == abs(max_entry_pivot_col))

    # If current pivot is not the max abs pivot value, swap rows:
    if (current_pivot != max_entry_pivot_col & pivot_row != m & pivot_col != n) {
      x = SwapRows(x, pivot_row, new_pivot_row)
    }

    # --------------------------------------------------------------------------
    # Eliminate rows
    # --------------------------------------------------------------------------
    # Define a counter to loop through rows below the pivot row
    start_elimination = pivot_row + 1

    # As long as the pivot row is not the last row, proceed with row elimination:
    if (pivot_row != m) {
      for (current_row in start_elimination:m) {
        x = EliminateRows(x, pivot_row, pivot_col, current_row)

        # Keep track of FLOPs per eliminated row
        FLOPs = FLOPs + (2*n + 1)
      }
    }

    # Continue while loop
    pivot_row = pivot_row + 1
    pivot_col = pivot_col + 1
  }

  # Generate output
  output = list(FLOPs = paste0("FLOPs: ", FLOPs), Echelon_Matrix = x)

  return(output)
}

# ==============================================================================
# Test case: 11 x 10
# ==============================================================================
# Generate 11 x 10 matrix with random uniform (0, 1) entries
set.seed(1)
test_matrix = matrix(data = runif((11*10), min = 0, max = 1), nrow = 11, ncol = 10)
test_matrix
SimpleEchelonElimination(test_matrix)

# ==============================================================================
# Test case: square
# ==============================================================================
square_matrix = matrix(data = runif((6*6), min = 0, max = 1), nrow = 6, ncol = 6)
square_matrix
SimpleEchelonElimination(square_matrix)

# ==============================================================================
# Test case: more columns than rows
# ==============================================================================
test_matrix_2 = matrix(data = runif((4*7), min = 0, max = 1), nrow = 4, ncol = 7)
test_matrix_2
SimpleEchelonElimination(test_matrix_2)
