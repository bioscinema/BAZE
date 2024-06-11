#' Detailed Taxonomic Tree Plot with Highlighted Nodes
#'
#' This function generates a detailed taxonomic tree plot with highlighted nodes and customized labels.
#' The function takes a phyloseq object and annotation data to create a tree plot, highlighting selected nodes,
#' and annotating them with specified labels.
#'
#' @param ps A phyloseq object containing the microbiome data.
#' @param anno.data A data frame containing annotation data with columns for nodes and colors.
#' @param alpha A numeric value specifying the transparency level for highlighting nodes. Default is 0.2.
#' @param anno.depth An integer specifying the depth of annotation. Default is 3.
#' @param anno.x An integer specifying the x-coordinate offset for annotations. Default is 10.
#' @param anno.y An integer specifying the y-coordinate offset for annotations. Default is 40.
#' @param layout A character string specifying the layout of the tree plot. Default is "circular".
#' 
#' @return A ggplot object representing the annotated taxonomic tree.
#' 
#' @examples
#' \dontrun{
#' ps <- phyloseq_object  # Replace with actual phyloseq object
#' anno.data <- data.frame(node = c("node1", "node2"), color = c("red", "blue"))  # Replace with actual annotation data
#' plot <- detailed_tax(ps, anno.data, alpha=0.2, anno.depth=3, layout="rectangular")
#' print(plot)
#' }
#' 
#' @importFrom dplyr arrange
#' @importFrom tidytree offspring
#' @importFrom ggtree ggtree geom_hilight geom_point2 geom_text2
#' @importFrom ggplot2 aes theme element_blank guides guide_legend
#' 
#' @export
detailed_tax <- function(ps, anno.data, alpha=0.2, anno.depth=3, anno.x=10, anno.y=40, layout="circular") {
  ## Generate taxonomy tree
  select_id <- get_selected_id(ps, anno.data)
  ps1 <- prune_taxa(select_id, ps)
  
  ## Remove duplicate taxa
  ps1 <- fix_duplicate_tax(ps1)
  
  ## Convert phyloseq to tree data
  tr1 <- phy_to_tax(ps1)
  
  plottree <- function(tree, size=0.5, layout=layout, shape=21, fill='white', color='black') {
    ggtree(tree, size=size, layout=layout)
  }
  gtree <- plottree(tr1)
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
  
  ## Add highlight
  for (i in 1:length(node_ids)) {
    n <- node_ids[i]
    color <- hilight.color[i]
    anno[n] <- color
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    offset <- get_offset(nodeClass)
    gtree <- gtree + geom_hilight(node=n, fill=color, alpha=alpha, extend=offset)
  }
  
  gtree$layers <- rev(gtree$layers)
  
  ## Add points and text labels to selected nodes with color
  for (i in 1:length(node_ids)) {
    n <- node_ids[i]
    color <- hilight.color[i]
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    if (nodeClass <= anno.depth) { ## Species and strains
      lab <- short.labs[1]
      short.labs <- short.labs[-1]
      if (is.null(short.labs.anno)) {
        short.labs.anno <- data.frame(lab=lab, annot=mapping$label, stringsAsFactors=FALSE)
      } else {
        short.labs.anno <- rbind(short.labs.anno, data.frame(lab=lab, annot=mapping$label, stringsAsFactors=FALSE))
      }
    } else {
      lab <- mapping$label
    }
    offset <- get_offset(nodeClass) - 0.4
    angle <- get_angle(n) + 90
    
    gtree <- gtree +
      geom_point2(data = data.frame(node = n, label = lab, x = mapping$x, y = mapping$y), 
                  aes(x = x, y = y), 
                  size = 8 + sqrt(nodeClass), 
                  color = color, 
                  fill = color, 
                  shape = 21, 
                  show.legend = FALSE) +
      geom_text2(data = data.frame(node = n, label = lab, x = mapping$x, y = mapping$y), 
                 aes(x = x, y = y, label = label), 
                 size = 5 + sqrt(nodeClass),
                 color="white",
                 vjust = 0.5, 
                 hjust = 0.5)
  }
  
  if (is.null(short.labs.anno)) { return(gtree) }
  
  ## Add short labels
  anno_shapes <- sapply(short.labs.anno$lab, utf8ToInt)
  gtree <- gtree + geom_point(data=short.labs.anno,
                              aes(x=0, y=0, shape=factor(annot)),
                              size=0, stroke=0) +
    guides(
      shape=guide_legend(override.aes=list(size=3, shape=anno_shapes))
    ) +
    theme(legend.position="right",
          legend.title=element_blank())
  
  return(gtree)
}
