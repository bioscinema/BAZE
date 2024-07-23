#' Visualize Taxonomic Tree with Annotations
#'
#' The `plottax` function generates a circular phylogenetic tree visualization with clades highlighted and labeled based on specified taxonomic annotations.
#'
#' @param tree A phylogenetic tree object compatible with `ggtree`, representing the evolutionary relationships among the taxa.
#' @param anno.data A data frame containing annotation information for the tree. It should include columns for `node` and `color`.
#' @param alpha A numeric value specifying the transparency level for the highlighted clades. Default is 0.2.
#' @param anno.depth A numeric value specifying the depth of annotation. Only nodes at or below this depth will be annotated with short labels. Default is 3.
#' @param anno.x A numeric value specifying the x-coordinate for annotation placement. Default is 10.
#' @param anno.y A numeric value specifying the y-coordinate for annotation placement. Default is 40.
#'
#' @details
#' The function proceeds through several steps to generate the final plot:
#' \itemize{
#'   \item Generates the initial ggtree plot.
#'   \item Highlights clades based on the provided annotation data.
#'   \item Adds labels to the clades based on their taxonomic level.
#'   \item Adds a legend for the annotations.
#' }
#'
#' @return Returns a ggtree plot object with highlighted clades and labels based on the specified taxonomic annotations.
#'
#' @importFrom ggtree ggtree geom_point geom_text2 geom_hilight
#' @importFrom dplyr filter arrange
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
#' p <- plottax(tree, anno.data, alpha = 0.2, anno.depth = 3, anno.x = 10, anno.y = 40)
#' print(p)
#' }
#'
#' @export
plottax <- function(tree, anno.data, alpha=0.2, anno.depth=3, anno.x=10, anno.y=40){
  # Generate taxonomy tree
  plottree <- function(tree, size=0.5, layout='circular', shape=21, fill='white', color='black'){
    ggtree(tree, size=size, layout=layout) +
      geom_point(aes(size=I(nodeSize)), shape=shape, fill=fill, color=color)
  }

  gtree <- plottree(tree)
  short.labs <- get_unique_id(length(unique(anno.data$node)))

  get_offset <- function(x) {(x*0.2 + 0.2)^2}

  get_angle <- function(node){
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

  for(i in 1:length(node_ids)){
    n <- node_ids[i]
    color <- hilight.color[i]
    anno[n] <- color
    offset <- get_offset(as.numeric(gtree$data[gtree$data$node == n, "nodeClass"]))
    gtree <- gtree + geom_hilight(node=n, fill=color, alpha=alpha, extend=offset)
  }

  gtree$layers <- rev(gtree$layers)
  gtree <- gtree + geom_point2(aes(size=I(nodeSize)), fill=anno, shape=21)

  short.labs.anno <- NULL
  for(i in 1:length(node_ids)){
    n <- node_ids[i]
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    if(nodeClass <= anno.depth){
      lab <- short.labs[1]
      short.labs <- short.labs[-1]
      if(is.null(short.labs.anno)){
        short.labs.anno <- data.frame(lab=lab, annot=mapping$label, stringsAsFactors=FALSE)
      } else {
        short.labs.anno <- rbind(short.labs.anno, data.frame(lab=lab, annot=mapping$label, stringsAsFactors=FALSE))
      }
    } else {
      lab <- mapping$label
    }
    offset <- get_offset(nodeClass) - 0.4
    angle <- get_angle(n) + 90
    gtree <- gtree + geom_text2(data=data.frame(node=n, label=lab, x=mapping$x, y=mapping$y),
                                aes(x=x, y=y, label=label), size=5 + sqrt(nodeClass),
                                angle=angle, vjust=-0.5, hjust=0.5)
  }

  if(is.null(short.labs.anno)){ return(gtree) }

  # Combine short labels and annotations for the legend
  short.labs.anno$legend <- paste(short.labs.anno$lab, short.labs.anno$annot, sep=": ")

  gtree + geom_point(data=short.labs.anno, aes(x=0, y=0, shape=factor(lab), fill=annot), size=0, stroke=0) +
    guides(shape=guide_legend(override.aes=list(size=3), label.position="right"), fill=FALSE) +
    theme(legend.position="right", legend.title=element_blank()) +
    scale_shape_manual(values=rep(21, length(short.labs.anno$lab)), labels=short.labs.anno$legend)
}
