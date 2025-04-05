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

  cat("Optimal Sampling Results:\n")
  cat(rep("=", 80), "\n", sep = "")

  # Pre-process to find paired elements
  processed_elements <- character(0)

  for (i in seq_along(x)) {
    element_name <- names(x)[i]
    element <- x[[i]]

    # Skip label element and already processed elements
    if (element_name == "label" || element_name %in% processed_elements) next

    # Check for w/w0 pair
    if (grepl("^w", element_name)) {
      counterpart_name <- ifelse(element_name == "w", "w0", "w")
      if (counterpart_name %in% names(x)) {
        # Print w and w0 as a table
        cat("Optimal approximate allocation:\n")
        print_weight_table(x[[element_name]], x[[counterpart_name]], x$label)
        processed_elements <- c(processed_elements, element_name, counterpart_name)
        cat(rep("-", 80), "\n", sep = "")
        next
      }
    }

    # Check for allocation (with or without allocation.real)
    if (grepl("^allocation", element_name)) {
      if (element_name == "allocation") {
        if ("allocation.real" %in% names(x)) {
          # Print both allocation and allocation.real as table
          cat("Optimal exact allocation:\n")
          print_allocation_table(x$allocation, x$allocation.real, x$label)
          processed_elements <- c(processed_elements, "allocation", "allocation.real")
        } else {
          # Print just allocation as single-row table
          cat("Optimal exact allocation:\n")
          print_single_allocation_table(x$allocation, x$label)
          processed_elements <- c(processed_elements, "allocation")
        }
        cat(rep("-", 80), "\n", sep = "")
        next
      }
    }

    # Default printing for non-paired elements
    cat(element_name, ":\n")

    if (is.matrix(element)) {
      print(format(element, justify = "right"), quote = FALSE)
    } else if (is.data.frame(element)) {
      print(element, row.names = FALSE)
    } else if (is.numeric(element)) {
      formatted <- format_numeric(element)
      cat(paste(formatted, collapse = ", "), "\n")
    } else if (is.logical(element)) {
      cat(paste(element, collapse = ", "), "\n")
    } else if (is.character(element)) {
      cat(paste0('"', element, '"', collapse = ", "), "\n")
    } else if (is.list(element)) {
      cat("Nested list with", length(element), "elements\n")
    } else {
      cat("Unsupported element type:", class(element), "\n")
    }

    cat(rep("-", 80), "\n", sep = "")
  }

  invisible(x)
}

# Helper function to print weight tables
print_weight_table <- function(w, w0, labels) {
  if (is.null(labels)) labels <- seq_along(w)

  # Create the table
  tbl <- rbind(
    w = format_numeric(w),
    w0 = format_numeric(w0)
  )

  colnames(tbl) <- labels
  print(tbl, quote = FALSE)
}

# Helper function to print allocation tables
print_allocation_table <- function(allocation, allocation_real, labels) {
  if (is.null(labels)) labels <- seq_along(allocation)

  # Create the table
  tbl <- rbind(
    allocation = format_numeric(allocation),
    allocation.real = format_numeric(allocation_real)
  )

  colnames(tbl) <- labels
  print(tbl, quote = FALSE)
}


# Helper function for single allocation row
print_single_allocation_table <- function(allocation, labels) {
  if (is.null(labels)) labels <- seq_along(allocation)

  # Create single-row table
  tbl <- rbind(
    allocation = format_numeric(allocation)
  )

  colnames(tbl) <- labels
  print(tbl, quote = FALSE)
}


# Helper function to format numeric values consistently
format_numeric <- function(x) {
  sapply(x, function(value) {
    if (value == round(value, digits = 0) && grepl("\\.0$", format(value, nsmall = 1))) {
      format(value, nsmall = 1)
    } else if (value == round(value)) {
      as.character(value)
    } else if (abs(value) < 1e-3 || abs(value) >= 1e+3) {
      formatC(value, format = "e", digits = 4)
    } else {
      sub("\\.?0+$", "", formatC(value, format = "f", digits = 4))
    }
  })
}
