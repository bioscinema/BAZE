##' @title fix_duplicate_tax
##'
##' @param physeq a phyloseq object
##' @author Chenghao Zhu, Chenhao Li, Guangchuang Yu
##' @export
##' @description fix the duplicatae taxonomy names of a phyloseq object

fix_duplicate_tax = function(physeq){
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package \"phyloseq\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  taxtab <- phyloseq::tax_table(physeq)
  for(i in 3:ncol(taxtab)){
    uniqs = unique(taxtab[,i])
    for(j in 1:length(uniqs)){
      if(is.na(uniqs[j])) next
      ind = which(taxtab[,i]== as.character(uniqs[j]))
      if(length(unique(taxtab[ind,i-1]))>1){
        taxtab[ind,i] = paste(taxtab[ind,i-1], taxtab[ind,i], sep="_")
      }
    }
  }
  phyloseq::tax_table(physeq) = taxtab
  return(physeq)
}

######################################################################
##' @title summarize_taxa
##'
##' @param physeq a phyloseq object
##' @param level the taxonomy level to summarize
##' @importFrom magrittr "%>%"
##' @importFrom reshape2 melt dcast
##' @import dplyr
##' @export
##' @description Summarize a phyloseq object on different taxonomy level

summarize_taxa = function(physeq, level, keep_full_tax = TRUE){
  # do some checking here
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package \"phyloseq\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  otutab = phyloseq::otu_table(physeq)
  taxtab = phyloseq::tax_table(physeq)
  
  if(keep_full_tax){
    taxonomy = apply(taxtab[,1:level], 1, function(x)
      paste(c("r__Root", x), collapse="|"))
  }else{
    taxonomy = taxtab[,level]
  }
  
  otutab %>%
    as.data.frame %>%
    mutate(taxonomy = taxonomy) %>%
    melt(id.var = "taxonomy",
         variable.name = "sample_id") %>%
    group_by(taxonomy, sample_id) %>%
    summarize(value=sum(value)) %>%
    dcast(taxonomy~sample_id)
}

################################################################################

#' Create a Phylogenetic Tree from Taxonomy Profile
#'
#' The `createtree` function constructs a phylogenetic tree object from a taxonomy profile. The tree includes node sizes based on relative abundance and classifies nodes into taxonomic levels.
#'
#' @param tax.profile Either a character string specifying the path to a taxonomy profile file or a data frame containing the taxonomy profile.
#' @param index An integer specifying the column index of the taxonomy string in the profile. Default is 1.
#' @param header A logical value indicating whether the taxonomy profile file has a header. Default is FALSE.
#' @param delimiter A character string specifying the delimiter used to separate taxonomic levels in the taxonomy string. Default is '\\|'.
#' @param node.size.scale A numeric value specifying the scaling factor for node sizes based on relative abundance. Default is 1.
#' @param node.size.offset A numeric value specifying the offset for node sizes. Default is 1.
#'
#' @details
#' The function reads a taxonomy profile, processes it to remove unclassified taxa, and splits the taxonomy strings into different levels. It then constructs a phylogenetic tree, assigns node sizes based on relative abundance, and classifies nodes into taxonomic levels.
#'
#' @return Returns a `treedata` object that includes the phylogenetic tree and associated data such as node sizes and classifications.
#'
#' @importFrom dplyr filter mutate select
#' @importFrom treeio treedata
#' @importFrom utils read.table
#' 
#'
#' @examples
#' \dontrun{
#' # Example usage
#' tax.profile <- "path/to/taxonomy_profile.txt"
#' tree <- createtree(tax.profile, index = 1, header = FALSE, delimiter = '\\|')
#' }
#'
#' @export
createtree <- function(tax.profile, index=1, header=FALSE, delimiter='\\|', node.size.scale=1, node.size.offset=1){
  if (is.character(tax.profile)) {
    taxtab <- read.table(tax.profile, sep='\t', stringsAsFactors=FALSE, header=header)
  }else{
    taxtab <- tax.profile
  }
  names(taxtab)[index] <- 'tax'
  names(taxtab)[-index] <- 'rel_abun'
  taxtab$tax <- as.character(taxtab$tax)
  taxtab <- taxtab %>% dplyr::filter(!grepl('unclassified|uncultured', taxtab[[index]])) # remove unclassified taxa
  tax_chars <- c('k', 'p', 'c', 'o', 'f', 'g', 's', 't')
  tax_split <- strsplit(taxtab$tax, delimiter)    ## split into different taxonomy levels
  child <- vapply(tax_split, tail, n=1, '')
  tax_class <- do.call(rbind, strsplit(child, '__'))[,1]
  parent <- vapply(tax_split, function(x) ifelse(length(x)>1, x[length(x)-1], 'root'), '')
  isTip <- !child %in% parent
  index <- c()
  index[isTip] <- 1:sum(isTip)
  index[!isTip] <- (sum(isTip)+1):length(isTip)
  ## tips comes first
  mapping <- data.frame(node=index, row.names=child, isTip, taxaAbun=taxtab$rel_abun)
  edges <- cbind(mapping[parent,]$node, mapping$node)
  edges <- edges[!is.na(edges[,1]),]
  
  a <- node.size.scale
  b <- node.size.offset
  mapping$nodeSize <- a*log(mapping$taxaAbun) + b
  mapping$nodeClass <- factor(tax_class, levels = rev(tax_chars))
  
  mapping <- mapping[order(mapping$node),]
  
  node.label <- rownames(mapping)[!mapping$isTip]
  phylo <- structure(list(edge = edges,
                          node.label = node.label,
                          tip.label = rownames(mapping[mapping$isTip,]),
                          edge.length=rep(1, nrow(edges)),
                          Nnode = length(node.label)
  ),
  class = "phylo")
  
  d <- mapping %>% dplyr::select_(~-isTip)
  treeio::treedata(phylo = phylo, data = dplyr::as_data_frame(d))
}
########################################################################

#' Convert Phyloseq Object to Phylogenetic Tree Data
#'
#' The `phy_to_tax` function converts a phyloseq object into a phylogenetic tree data structure suitable for use with various phylogenetic analysis and visualization tools.
#'
#' @param physeq A phyloseq object containing the microbial community data.
#' @param use_abundance A logical value indicating whether to use abundance data. Default is TRUE.
#' @param node.size.scale A numeric value specifying the scaling factor for node sizes based on relative abundance. Default is 1.
#' @param node.size.offset A numeric value specifying the offset for node sizes. Default is 1.
#'
#' @details
#' The function processes the phyloseq object to ensure that taxa are rows, fixes any duplicate taxonomy, and summarizes the taxa at each level. It then constructs a tree data structure using the `createtree` function.
#'
#' @return Returns a `treedata` object that includes the phylogenetic tree and associated data such as node sizes and classifications.
#'
#' @importFrom phyloseq phyloseq tax_table otu_table taxa_are_rows
#' @importFrom dplyr filter mutate as_data_frame
#' @importFrom treeio treedata
#'
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(phyloseq)
#' library(ggtree)
#' 
#' # Example usage
#' data(GlobalPatterns)
#' tree_data <- phy_to_tax(GlobalPatterns)
#' }
#'
#' @export
phy_to_tax <- function(physeq,
                          use_abundance = TRUE,
                          node.size.scale = 1,
                          node.size.offset = 1){
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package \"phyloseq\" is needed for this function. Please install it.",
         call. = FALSE)
  }
  
  fix_duplicate_tax <- function(physeq){
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("Package \"phyloseq\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    taxtab <- phyloseq::tax_table(physeq)
    for(i in 3:ncol(taxtab)){
      uniqs = unique(taxtab[,i])
      for(j in 1:length(uniqs)){
        if(is.na(uniqs[j])) next
        ind = which(taxtab[,i]== as.character(uniqs[j]))
        if(length(unique(taxtab[ind,i-1]))>1){
          taxtab[ind,i] = paste(taxtab[ind,i-1], taxtab[ind,i], sep="_")
        }
      }
    }
    phyloseq::tax_table(physeq) = taxtab
    return(physeq)
  }
  
  # convert taxa_are_rows to TRUE if not
  if(!phyloseq::taxa_are_rows(physeq)){
    otu = phyloseq::otu_table(physeq)
    otu = phyloseq::otu_table(t(otu@.Data), taxa_are_rows = TRUE)
    phyloseq::otu_table(physeq) = otu
  }
  
  taxtab <- tryCatch(
    phyloseq::tax_table(physeq),
    error = function(e){
      stop("The tax_table is required to draw the cladogram")
    }
  )
  if (use_abundance) {
    otutab <- phyloseq::otu_table(physeq)
  } else {
    row.names = rownames(phyloseq::tax_table(physeq))
    otutab = matrix(rep(0, length(row.names)), ncol = 1)
    rownames(otutab) <- row.names
  }
  # Sometimes the upper level taxonomy is NA, for example:
  # r__Root|k__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Rickettsiales|NA|g__CandidatusPelagibacter
  # Remove all the labels after NA, because they not trustable.
  for(i in 2:ncol(taxtab)){
    na.idx.upper.taxa <- is.na(taxtab[,i])
    na.idx.current.taxa <- is.na(taxtab[, i - 1])
    idx.to.update <- (!na.idx.upper.taxa) & (na.idx.current.taxa)
    if (any(idx.to.update)) {
      taxtab[idx.to.update, i] = NA
    }
  }
  
  # add taxonomy level if not already have
  if(!grepl("^k__", taxtab[1,1])){
    for(i in 1:ncol(taxtab)){
      tax_level = tolower(strsplit(colnames(taxtab)[i],'')[[1]][1])
      taxtab[,i] = paste0(tax_level, "__", taxtab[,i])
    }
  }
  
  # summarize taxa
  if(use_abundance){
    otutab_2 = otutab %>%
      rowSums %>%
      data.frame()
    names(otutab_2) = "value"
  }else{
    otutab_2 = data.frame(
      row.names = rownames(otutab),
      value = rep(0, nrow(otutab))
    )
  }
  physeq_2 = phyloseq::phyloseq(
    otu_table(otutab_2, taxa_are_rows = TRUE),
    taxtab
  )
  
  treetable = data.frame(taxonomy = "r__Root", value = 0)
  if(use_abundance){
    treetable$value = 100
  }
  
  for(i in 1:ncol(taxtab)){
    summarized = summarize_taxa(physeq_2, level=i)
    summarized = tidytree::filter(summarized, !grepl("NA$", summarized$taxonomy))
    if(use_abundance){
      summarized = mutate(
        summarized,
        value = value/sum(value) * 100
      )
    }
    treetable = rbind(treetable, summarized)
  }
  if(!use_abundance) treetable$value = treetable$value + 5
  
  createtree(treetable,
                    index = 1,
                    header = FALSE,
                    delimiter = "\\|",
                    node.size.scale = node.size.scale,
                    node.size.offset = node.size.offset)
}


##########################################################################
#' Get Selected Taxonomic IDs
#'
#' This function filters the taxonomic table of a phyloseq object to retain only the taxa that match
#' the specified nodes in the annotation data.
#'
#' @param ps A phyloseq object containing the microbiome data.
#' @param anno.data A data frame containing annotation data with columns for nodes.
#' 
#' @return A vector of selected taxonomic IDs.
#' 
#' @examples
#' \dontrun{
#' ps <- phyloseq_object  # Replace with actual phyloseq object
#' anno.data <- data.frame(node = c("g__Bacteroides", "f__Lachnospiraceae"), color = c("red", "blue"))  # Replace with actual annotation data
#' selected_ids <- get_selected_id(ps, anno.data)
#' print(selected_ids)
#' }
#' 
#' @importFrom dplyr mutate filter across
#' @export
get_selected_id <- function(ps, anno.data) {
  tax_tab <- as.data.frame(tax_table(ps))
  
  # Add the appropriate prefixes to each taxonomic level
  tax_tab_prefixed <- tax_tab %>%
    mutate(across(everything(), ~ paste0(tolower(substr(cur_column(), 1, 1)), "__", .)))
  
  # Filter the tax table to keep only the selected features
  filtered_tax_tab <- tax_tab_prefixed %>%
    filter(apply(., 1, function(row) any(row %in% anno.data$node)))
  
  selected_id <- rownames(filtered_tax_tab)
  return(selected_id)
}
