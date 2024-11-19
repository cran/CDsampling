#' Print Method for matrix_list Objects
#'
#' @description Custom print method for objects of class `matrix_list`.
#' @param x An object of class `matrix_list`.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x` (the input object).
#' @export
#' @method print matrix_list
#'
print.matrix_list <- function(x, ...) {
  cat("\nF_x:\n")
  cat("Dimensions:", dim(x$F_x)[1], "x", dim(x$F_x)[2], "\n")
  separator1 <- paste(rep("-", ncol(x$F_x) * 5), collapse = "")
  cat(separator1, "\n")
  print(format(x$F_x, justify = "right"), quote = FALSE)
  cat(separator1, "\n")

  cat("\nU_x:\n")
  cat("Dimensions:", dim(x$U_x)[1], "x", dim(x$U_x)[2], "\n")
  separator2 <- paste(rep("-", ncol(x$U_x) * 5), collapse = "")
  cat(separator2, "\n")
  print(format(x$U_x, justify = "right"), quote = FALSE)
  cat(separator2, "\n")

  invisible(x)
}
