#' Create Taxonomy Annotation File
#'
#' This function generates an annotation file for a taxonomy tree based on a merged table of taxonomic data. It assigns colors to taxa based on the effect sign and writes various annotations to the specified file.
#'
#' @param merged_table A data frame containing taxonomic data. It must have columns `taxa` and `effect_sign`.
#' @param annotation_file A character string specifying the path to the output annotation file.
#' @details
#' The function first assigns colors to the taxa based on their `effect_sign` values: "red" for positive effects and "blue" for negative effects. It then writes general annotation settings to the specified annotation file. For each taxon in the merged table, if the taxon is not "unknown" and does not contain "uncultured", the function writes detailed annotations to the file, including marker size, color, rotation, background color, and font size.
#' @examples
#' # Example data frame
#' merged_table <- data.frame(
#'   taxa = c("Genus1", "uncultured_bacteria", "Genus2", "unknown", "uncultured_fungus", "Family1"),
#'   effect_sign = c("positive", "negative", "positive", "negative", "negative", "positive")
#' )
#' # Define the annotation file path
#' annotation_file <- "taxonomy_annotations.txt"
#' # Call the function
#' create_tax_merged(merged_table, annotation_file)
#' @export
create_tax_merged <- function(merged_table, annotation_file) {
  ## create color for different signs
  merged_table$color <- ifelse(merged_table$effect_sign == "positive", "red", "blue")
  
  file_conn <- file(annotation_file, open = "wt")
  ## write general annotation first
  writeLines("title\tTaxonomy tree with selection result", file_conn)
  writeLines("title_font_size\t25", file_conn)
  writeLines("start_rotation\t90", file_conn)
  # writeLines("clade_separation\t0.35", file_conn)
  writeLines("class_legend_font_size\t10", file_conn)
  writeLines("annotation_legend_font_size\t15", file_conn)
  writeLines("annotation_font_size\t10", file_conn)
  writeLines("clade_marker_color\twhite", file_conn)
  writeLines("branch_color_from_ancestor\t0", file_conn)
  # writeLines("annotation_background_separation\t-0.24", file_conn)
  # writeLines("annotation_background_offset\t-0.21", file_conn)
  # writeLines("annotation_background_width\t0.03", file_conn)
  writeLines("branch_thickness\t1.5", file_conn)
  writeLines("branch_bracket_depth\t0.5", file_conn)
  
  for (i in 1:nrow(merged_table)) {
    otu <- merged_table$taxa[i]
    color <- merged_table$color[i]
    
    if (otu == "unknown" || 
        grepl("uncultured", otu, ignore.case = TRUE) || 
        grepl("unclassified", otu, ignore.case = TRUE) || 
        grepl("unclass", otu, ignore.case = TRUE)) {
      next
    }
    writeLines(paste(otu, "clade_marker_size", "300", sep = "\t"), file_conn)
    writeLines(paste(otu, "clade_marker_color", color, sep = "\t"), file_conn)
    writeLines(paste(otu, "annotation_rotation", "270", sep = "\t"), file_conn)
    writeLines(paste(otu, "annotation_background_color", "white", sep = "\t"), file_conn)
    writeLines(paste(otu, "annotation", otu, sep = "\t"), file_conn)
    writeLines(paste(otu, "annotation_font_size", "18", sep = "\t"), file_conn)
  }
  
  close(file_conn)
}
