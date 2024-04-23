#' Create Taxonomic Annotations for Visualization
#'
#' This function generates annotation files from a phyloseq object based on specified taxonomic levels.
#' It validates the taxonomic structure, computes relative abundances, and assigns colors to different
#' taxonomic levels for visualization purposes. The annotations are saved to a file which can be used in
#' various tree visualization tools to enhance the presentation of phylogenetic data.
#'
#' @param ps A `phyloseq` object containing microbiome or other taxonomic data.
#' @param annotation_file Character, path to the file where annotations will be saved.
#' @param level Character, the taxonomic level at which annotations are to be generated.
#'        Default is "Genus". Supported levels depend on the taxonomy table in the `phyloseq` object.
#'
#' @return None; this function writes directly to a file specified by `annotation_file`.
#'
#' @details
#' The function first checks if taxa are correctly set as rows in the `phyloseq` object. It then extracts
#' the OTU table and taxonomy table, identifies the unique taxa at the specified level, and assigns a
#' unique color to each. It calculates the relative abundance of each taxon, scales these values, and
#' writes detailed annotation data into the specified file. This file includes settings for colors, sizes,
#' and other visual attributes for each taxon, facilitating detailed customization of phylogenetic tree visualizations.
#'
#' @examples
#' \dontrun{
#'   library(phyloseq)
#'   data(GlobalPatterns)
#'   create_tax_annot(GlobalPatterns, "path/to/your/annotation_file.txt", level="Genus")
#' }
#'
#' @importFrom phyloseq otu_table tax_table taxa_are_rows
#' @importFrom stats aggregate
#' @importFrom randomcoloR distinctColorPalette
#' @export
create_tax_annot <- function(ps, annotation_file, level="Genus"){
  ## Check if taxa are set as rows in the phyloseq object
  if(!taxa_are_rows(ps)){
    stop("please check your phyloseq subject and make sure your taxa are rows")
  }

  ## Extract otu table and taxonomy table
  myotu <- as.data.frame(otu_table(ps))
  mytax <- as.data.frame(tax_table(ps))
  otu_all <- row.names(myotu)
  expected_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  actual_ranks <- colnames(mytax)
  missing_ranks <- expected_ranks[!expected_ranks %in% actual_ranks]
  if (length(missing_ranks) > 0) {
    stop(paste("The following taxonomic rank names are missing or do not match exactly:", paste(missing_ranks, collapse=", "), ". Please check your column names."))
  }
  unique_phyla <- unique(mytax$Phylum)
  num_colors <- length(unique_phyla)
  phylum_colors <- distinctColorPalette(num_colors)
  color_list <- as.data.frame(cbind(phylum_colors, unique_phyla))

  ##calculae relative abundance
  merged_table <- merge(mytax, myotu, by="row.names", all = TRUE)
  colnames(merged_table[1]) <- c("OTU")
  myotu$expected_level <- mytax[,colnames(mytax)==level]
  level_abundance <- aggregate(. ~ expected_level, data = myotu, sum)
  level_abundance$total <- rowSums(level_abundance[,-1])/sum(level_abundance[,-1])
  level_abundance$size <- level_abundance$total * 100000

  file_conn <- file(annotation_file, open = "wt")
  writeLines("title\tTaxonomy tree", file_conn)
  writeLines("title_font_size\t25", file_conn)
  writeLines("annotation_legend_font_size\t15", file_conn)
  writeLines("annotation_font_size\t15", file_conn)
  writeLines("branch_thickness\t2", file_conn)
  writeLines("branch_bracket_depth\t0.5", file_conn)
  writeLines("branch_color_from_ancestor\t0", file_conn)

  for (otu in otu_all) {
    Species <- as.character(mytax[otu,]$Species)
    Genus <- as.character(mytax[otu,]$Genus)
    Class <- as.character(mytax[otu,]$Class)
    Order <- as.character(mytax[otu,]$Order)
    Family <- as.character(mytax[otu,]$Family)
    Kingdom <- as.character(mytax[otu,]$Kingdom)
    Phylum <- as.character(mytax[otu,]$Phylum)
    level_tax <- as.character(mytax[otu,colnames(mytax)==level])

    first_letter_class <- toupper(substr(Class, 1, 1))
    first_letter_phylum <- toupper(substr(Phylum,1,1))
    if (Species == "unknown" | Class == "unknown" | Order == "unknown" | Family == "unknown" |
        Kingdom == "unknown" | Phylum == "unknown" | Genus == "unknown" | level_tax == "unknown") {
      next
    }

    color <- color_list$phylum_colors[color_list$unique_phyla==Phylum]
    size <- ifelse(!is.na(level_abundance$size[level_abundance$expected_level == level_tax]),
                   level_abundance$size[level_abundance$expected_level == level_tax],
                   0)
    # writeLines(paste(paste(kingdom, phylum, class, order, family,genus,species, sep = "."),"clade_marker_size","40", sep = "\t"), file_conn)
    # writeLines(paste(paste(kingdom, phylum, class, order, family,genus,species, sep = "."), "clade_marker_color", color, sep = "\t"), file_conn)
    writeLines(paste(level_tax, "ring_height","1",size/10000, sep = "\t"), file_conn)
    writeLines(paste(level_tax,"ring_color", "1","purple", sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, sep = "."), "clade_marker_size", "100",sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, sep = "."), "clade_marker_color", color,sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, sep = "."), "annotation_background_color", color,sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, sep = "."), "annotation",Phylum, sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, sep = "."), "annotation_rotation","90", sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, sep = "."),"annotation_font_size","10", sep = "\t"), file_conn)
    # writeLines(paste(paste(Kingdom, Phylum, Class,sep = "."),"annotation", paste0(first_letter_class,":",Class), sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, Class,sep = "."),"annotation_background_color", color, sep = "\t"), file_conn)
    # writeLines(paste(class, "annotation_rotation","90", sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, Class,sep = "."),"clade_marker_size", "90", sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, Class,sep = "."),"clade_marker_color", color, sep = "\t"), file_conn)
    # writeLines(paste(genus,"annotation", genus, sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, Class, Order, Family,Genus, sep = "."),"annotation_background_color", color, sep = "\t"), file_conn)
    # writeLines(paste(genus, "annotation_rotation","90", sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, Class, Order, Family,Genus, sep = "."),"clade_marker_size", "60", sep = "\t"), file_conn)
    writeLines(paste(paste(Kingdom, Phylum, Class, Order, Family,Genus, sep = "."),"clade_marker_color", color, sep = "\t"), file_conn)
    writeLines(paste(Order,"clade_marker_color",color,sep="\t"),file_conn)
    writeLines(paste(Order, "clade_marker_size","80",sep = "\t"),file_conn)
    writeLines(paste(Family,"clade_marker_color",color,sep = "\t"),file_conn)
    writeLines(paste(Family,"clade_marker_size","70",sep = "\t"),file_conn)
    writeLines(paste(Species,"clade_marker_size","50",sep = "\t"),file_conn)
    writeLines(paste(Species,"clade_marker_color",color,sep = "\t"),file_conn)
    }
  close(file_conn)
}
