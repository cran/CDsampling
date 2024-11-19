#' Print Method for matrix_output Objects
#'
#' @description Custom print method for objects of class `matrix_output`.
#' @param x An object of class `matrix_output`.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x` (the input object).
#' @export
#' @method print matrix_output
print.matrix_output <- function(x, ...) {

  cat("Dimensions:", dim(x)[1], "x", dim(x)[2], "\n")
  cat("Matrix:\n")

  # Add a separator for readability
  separator <- paste(rep("-", ncol(x) * 5), collapse = "")
  cat(separator, "\n")

  # Format the matrix with aligned columns
  print(format(x, justify = "right"), quote = FALSE)
  cat(separator, "\n")

  invisible(x)
}
