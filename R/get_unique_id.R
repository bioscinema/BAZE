#' Generate Unique Identifiers
#'
#' This function generates a specified number of unique identifiers
#' by combining letters from the English alphabet. It starts with single letters
#' and expands to combinations of multiple letters as needed.
#'
#' @param n Integer. The number of unique identifiers required.
#' @param depth Integer. The current depth of the identifier combinations (default is 1).
#' @return A character vector of unique identifiers.
#' @examples
#' get_unique_id(30)
#' @export
get_unique_id <- function(n, depth = 1) {
  # Create a list of 'depth' elements, each containing the alphabet letters
  args <- lapply(seq_len(depth), FUN = function(x) letters)
  
  # Generate all combinations of the letters up to the specified depth
  x <- do.call(expand.grid, args = list(args, stringsAsFactors = FALSE))
  
  # Reverse the column order for proper concatenation
  x <- x[, rev(names(x)), drop = FALSE]
  
  # Concatenate the letter combinations into strings
  x <- do.call(paste0, x)
  
  # Check if the required number of identifiers is generated
  if (n <= length(x)) {
    return(x[seq_len(n)])
  }
  
  # If more identifiers are needed, recursively call the function with increased depth
  return(c(x, get_unique_id(n - length(x), depth = depth + 1)))
}