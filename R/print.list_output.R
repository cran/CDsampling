#' Print Method for list_output Objects
#'
#' @description Custom print method for objects of class `list_output`.
#' @param x An object of class `list_output`.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x` (the input object).
#' @export
#' @method print list_output
#'
#'

print.list_output <- function(x, ...) {
  if (!is.list(x)) {
    stop("The object is not a list.")
  }

  cat("List Object:\n")
  cat("Number of elements:", length(x), "\n")

  for (i in seq_along(x)) {
    element <- x[[i]]
    element_name <- names(x)[i]

    cat("\nElement", i, if (!is.null(element_name)) paste0(" ('", element_name, "')"), ":\n")
    cat(rep("-", 60), "\n", sep = "")

    if (is.matrix(element)) {
      # Pretty print for matrices
      print(format(element, justify = "right"), quote = FALSE)
    } else if (is.data.frame(element)) {
      # Pretty print for data frames
      print(element, row.names = FALSE)
    } else if (is.atomic(element)) {
      # Print atomic vectors (e.g., numeric, character)
      cat(paste(element, collapse = ", "), "\n")
    } else if (is.list(element)) {
      # Indicate nested lists
      cat("Nested list with", length(element), "elements\n")
    } else {
      # Fallback for unknown types
      cat("Unsupported element type:", class(element), "\n")
    }

    cat(rep("-", 60), "\n", sep = "")
  }

  invisible(x)
}
