#' Visualize Phylogenetic Tree with Taxonomic Annotations and Relative Abundance Bars
#'
#' The `tax_rel` function generates a circular phylogenetic tree visualization with clades highlighted, annotated based on specified taxonomic levels, and includes bars representing relative abundance outside the tree.
#'
#' @param tree A phylogenetic tree object compatible with `ggtree`, representing the evolutionary relationships among the taxa.
#' @param anno.data A data frame containing annotation information for the tree. It should include columns for `node` and `color`.
#' @param alpha A numeric value specifying the transparency level for the highlighted clades. Default is 0.2.
#' @param anno.depth A numeric value specifying the depth of annotation. Only nodes at or below this depth will be annotated with short labels. Default is 3.
#' @param anno.x A numeric value specifying the x-coordinate for annotation placement. Default is 10.
#' @param anno.y A numeric value specifying the y-coordinate for annotation placement. Default is 40.
#' @param scale_size A numeric value to scale the size of the relative abundance bars and tree plot. Default is 0.5.
#'
#' @details
#' The function proceeds through several steps to generate the final plot:
#' \itemize{
#'   \item Generates the initial ggtree plot.
#'   \item Highlights clades based on the provided annotation data.
#'   \item Propagates phylum information from internal nodes to tips.
#'   \item Adds labels to the clades based on their taxonomic level.
#'   \item Adds bars representing relative abundance outside the tree.
#'   \item Adds a legend for the annotations.
#' }
#'
#' @return Returns a ggtree plot object with highlighted clades, labels, and relative abundance bars based on the specified taxonomic annotations.
#'
#' @importFrom ggtree ggtree geom_point geom_text2 geom_hilight
#' @importFrom ggplot2 geom_segment
#' @importFrom dplyr filter arrange mutate
#' @importFrom scales hue_pal
#'
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(ggtree)
#' library(dplyr)
#'
#' # Example usage
#' tree <- read.tree("path/to/tree_file") # Load your tree
#' anno.data <- data.frame(node = c("node1", "node2"), color = c("red", "blue"))
#' p <- tax_rel(tree, anno.data, alpha = 0.2, anno.depth = 3, anno.x = 10, anno.y = 40, scale_size = 0.5)
#' print(p)
#' }
#'
#' @export
tax_rel <- function(tree, anno.data, alpha=0.2, anno.depth=3, anno.x=10, anno.y=40, scale_size=0.5) {
  ## Generate taxonomy tree
  plottree <- function(tree, size=0.5, layout='circular', shape=21, fill='white', color='black') {
    ggtree(tree, size=size, layout = layout) +
      geom_point(aes(size=I(nodeSize)), shape=shape, fill=fill, color=color, show.legend = FALSE)
  }

  gtree <- plottree(tree)
  short.labs <- get_unique_id(length(unique(anno.data$node)))

  get_offset <- function(x) { (x * 0.2 + 0.2)^2 }

  get_angle <- function(node) {
    data <- gtree$data
    sp <- tidytree::offspring(data, node)$node
    sp2 <- c(sp, node)
    sp.df <- data[match(sp2, data$node),]
    mean(range(sp.df$angle))
  }

  anno.data <- arrange(anno.data, node)
  hilight.color <- anno.data$color
  node_list <- anno.data$node
  node_ids <- (gtree$data %>% filter(label %in% node_list) %>% arrange(label))$node
  anno <- rep('white', nrow(gtree$data))

  ## Extract phylum information from the labels
  gtree$data$phylum <- ifelse(grepl("p__", gtree$data$label), sub(".*(p__[^;]+).*", "\\1", gtree$data$label), NA)
  ## Propagate phylum information from internal nodes to tips
  internal_nodes <- gtree$data %>% filter(!isTip & !is.na(phylum)) %>% pull(node)
  for (node in internal_nodes) {
    descendant_nodes <- offspring(gtree$data, node)$node
    tip_nodes <- descendant_nodes[descendant_nodes %in% gtree$data$node[gtree$data$isTip]]
    gtree$data$phylum[gtree$data$node %in% tip_nodes] <- gtree$data$phylum[gtree$data$node == node]
  }

  ## Add highlight ... duplicated code
  for (i in 1:length(node_ids)) {
    n <- node_ids[i]
    color <- hilight.color[i]
    anno[n] <- color
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    offset <- get_offset(nodeClass)
    gtree <- gtree + geom_hilight(node = n, fill = color, alpha = alpha, extend = offset)
  }

  gtree$layers <- rev(gtree$layers)

  ## Set nodeSize to NA for all nodes initially
  gtree$data$nodeSize <- NA

  ## Set nodeSize for selected nodes
  selected_node_size <- 2  # Change this value as needed
  gtree$data$nodeSize[gtree$data$node %in% node_ids] <- selected_node_size

  ## Add points to selected nodes only
  gtree <- gtree + geom_point2(aes(size = I(nodeSize)), fill = anno, shape = 21, show.legend = FALSE)

  ## Add labels
  short.labs.anno <- NULL
  for (i in 1:length(node_ids)) {
    n <- node_ids[i]
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    if (nodeClass <= anno.depth) { ## Species and strains
      lab <- short.labs[1]
      short.labs <- short.labs[-1]
      if (is.null(short.labs.anno)) {
        short.labs.anno <- data.frame(node = n, lab = lab, annot = mapping$label, stringsAsFactors = FALSE)
      } else {
        short.labs.anno <- rbind(short.labs.anno, data.frame(node = n, lab = lab, annot = mapping$label, stringsAsFactors = FALSE))
      }
    } else {
      lab <- mapping$label
    }
    offset <- get_offset(nodeClass) - 0.4
    angle <- get_angle(n) + 90
  }

  if (is.null(short.labs.anno)) { return(gtree) }

  ## Combine short labels and annotations for the legend
  short.labs.anno$legend <- paste(short.labs.anno$lab, short.labs.anno$annot, sep=": ")

  gtree <- gtree + geom_point(data=short.labs.anno, aes(x=0, y=0, shape=factor(lab), fill=annot), size=0, stroke=0) +
    guides(shape=guide_legend(override.aes=list(size=3), label.position="right"), fill=FALSE) +
    theme(legend.position="right", legend.title=element_blank()) +
    scale_shape_manual(values=rep(21, length(short.labs.anno$lab)), labels=short.labs.anno$legend)

  if (!is.null(gtree$data$taxaAbun)) {
    tip_nodes <- gtree$data %>% filter(isTip) %>% filter(!is.na(taxaAbun))
    relative_abundance <- tip_nodes %>%
      mutate(xstart = max(gtree$data$x) + 1,  # Positioning the bars outside the circle
             xend = xstart + taxaAbun * scale_size)  # Scale the relative abundance for better visualization

    gtree <- gtree +
      geom_segment(data = relative_abundance,
                   aes(x = xstart, xend = xend, y = y, yend = y, color = phylum),
                   size = 2, show.legend = FALSE) +
      guides(color="none")
  }

  return(gtree)
}

# Usage example (replace with actual tree and anno.data)
# p2 <- plottax(tr1, anno.data = effect, anno.depth = 6)
# print(p2)
