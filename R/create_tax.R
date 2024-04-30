#' Export Cleaned Taxonomic Paths to a File
#'
#' This function processes taxonomic data from a `phyloseq` object, cleans it by removing unwanted
#' characters and filtering out non-informative entries, and writes the resulting taxonomic paths
#' to a specified file. Each path is a concatenation of the taxonomic levels for each organism,
#' providing a clear hierarchical view of the taxonomy.
#'
#' @param ps A `phyloseq` object containing taxonomic data.
#' @param tree_file Character, the filepath where the cleaned taxonomic paths will be written.
#'
#' @return None; the function writes the cleaned taxonomic paths directly to the specified file.
#'
#' @details
#' The function verifies that the taxonomic data are structured with taxa as rows. It then converts
#' the taxonomic table to a dataframe and cleanses it by removing periods and ignoring entries labeled
#' as "unknown" or NA. The cleaned taxonomic names are then concatenated to form paths, which are
#' written to the provided file, each path on a new line. This processing is essential for ensuring
#' the clarity and usability of taxonomic data in downstream analyses or visualizations.
#'
#' @examples
#' \dontrun{
#'   library(phyloseq)
#'   data(GlobalPatterns)
#'   create_tax(GlobalPatterns, "path/to/your/clean_taxonomic_paths.txt")
#' }
#'
#' @importFrom phyloseq tax_table taxa_are_rows
#' @export
create_tax <- function(ps, tree_file) {
  # Check if the phyloseq object has a taxonomic table
  if (!taxa_are_rows(ps)) {
    stop("The provided phyloseq object does not contain a taxonomic table.")
  }

  # Convert the taxonomic table to a data frame
  mytax <- as.data.frame(tax_table(ps))

  # Remove periods from taxonomic names
  mytax[] <- lapply(mytax, function(col) gsub("\\.", "", col))

  # Create taxonomic paths by concatenating taxonomic levels, omitting "unknown" and NAs
  tax_paths <- apply(mytax, 1, function(row) {
    # Convert all entries to lower case and remove undesired taxa
    # Filter out 'unknown', 'uncultured', and variations like 'uncultured_xxx'
    valid_taxa <- row[!grepl("^unknown$|^uncultured", tolower(row), perl = TRUE) & !is.na(row)]

    # Concatenate the remaining taxonomic information into a single string
    paste(na.omit(valid_taxa), collapse = ".")
  })


  # Write the taxonomic paths to the specified file
  writeLines(tax_paths, con = tree_file)
}
