#' Visualize Phylogenetic Tree with Taxonomic Clade Highlights
#'
#' The `taxview` function generates a phylogenetic tree visualization with clades highlighted and labeled based on a specified taxonomic level from a phyloseq object.
#' This function facilitates the interpretation of microbial community data by incorporating taxonomic information directly into the tree plot.
#'
#' @param ps A phyloseq object containing the microbial community data, including the OTU table and taxonomy table.
#' @param tree A phylogenetic tree object compatible with `ggtree`, representing the evolutionary relationships among the taxa in the phyloseq object.
#' @param branch_thickness A numeric value specifying the thickness of the branches in the phylogenetic tree. Default is 0.5.
#' @param layout A character string specifying the layout of the phylogenetic tree. Default is 'circular'. Other options include 'rectangular', 'slanted', etc.
#' @param level A character string specifying the taxonomic level to highlight in the tree (e.g., "Phylum", "Genus"). Default is "Phylum".
#'
#' @details
#' The function proceeds through several steps to generate the final plot:
#' \itemize{
#'   \item Extracts the taxonomy table from the phyloseq object.
#'   \item Generates a color palette for the unique taxonomic groups at the specified level.
#'   \item Maps the tips of the phylogenetic tree to their corresponding taxonomic groups.
#'   \item Creates the initial ggtree plot with the specified branch thickness and layout.
#'   \item Highlights and labels clades for each unique taxonomic group:
#'     \itemize{
#'       \item Identifies the internal node label in the tree corresponding to each taxonomic group.
#'       \item Highlights the clade associated with that node using the assigned color.
#'       \item Adds clade labels outside the highlighted clades.
#'     }
#'   \item Adds a legend mapping the taxonomic groups to their colors.
#' }
#'
#' @return Returns a ggtree plot object with highlighted clades and labels based on the specified taxonomic level.
#'
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(phyloseq)
#' library(ggtree)
#' library(scales)
#'
#' # Assuming `ps1` is your phyloseq object and `tree` is the corresponding phylogenetic tree
#' p <- taxview(ps1, tree, branch_thickness = 0.5, layout = 'circular', level = "Phylum")
#'
#' # Display the plot
#' print(p)
#' }
#'
#' @import phyloseq
#' @import ggtree
#' @import scales
#' 
#' @export
taxview <- function(ps, tree, branch_thickness=0.5, layout='circular', level="Phylum") {
  # Extract taxonomy table
  tax <- as.data.frame(ps@tax_table)
  
  # Generate a color palette for the taxonomic levels
  taxon_colors <- hue_pal()(length(unique(tax[[level]])))
  names(taxon_colors) <- unique(tax[[level]])
  
  # Create a mapping from tips to their taxonomy level
  tip_labels <- data.frame(label = rownames(tax), taxon = tax[[level]])
  
  # Create the ggtree plot
  p <- ggtree(tree, size=branch_thickness, layout=layout)
  
  # Highlight each taxon clade based on internal node labels and add clade labels
  for (taxon in unique(tip_labels$taxon)) {
    node_label <- paste0(tolower(substr(level, 1, 1)), "__", taxon)
    if (node_label %in% tree@phylo[["node.label"]]) {
      node_index <- which(tree@phylo[["node.label"]] == node_label) + length(tree@phylo[["tip.label"]])
      p <- p + geom_hilight(node = node_index, fill = taxon_colors[taxon], alpha = 0.3, show.legend = TRUE)
      
      # Add clade labels outside the highlighted clades
      p <- p + geom_cladelabel(node = node_index, label = taxon,
                               fontsize = 3, offset = 1, barsize = 0, hjust = 0.5, angle = 0)
    }
  }
  
  # Add a dummy dataframe to create the legend
  legend_df <- data.frame(taxon = names(taxon_colors), fill = taxon_colors)
  
  # Add the legend manually
  p <- p + geom_point(data = legend_df, aes(x = 0, y = 0, fill = taxon), size = 0) +
    scale_fill_manual(values = taxon_colors, name = level) +
    guides(fill = guide_legend(override.aes = list(size = 5, shape = 21))) +
    theme(legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  # Return the plot
  return(p)
}