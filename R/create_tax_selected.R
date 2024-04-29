#' Generate Detailed Taxonomic Annotations for Selected Taxa
#'
#' This function processes a `phyloseq` object to create annotations for taxa that are notably represented
#' in the results of variable selection model. It extracts taxonomic and abundance data, evaluates the
#' significance of each taxa's representation in the results, and writes annotations with specific visual
#' attributes to a file.
#'
#' @param ps A `phyloseq` object containing taxonomic data.
#' @param nburnin Integer, number of burn-in iterations to disregard in the analysis.
#' @param niter Integer, number of iterations considered for analysis after burn-in.
#' @param result List, containing statistical results including 'betahat' and 'gammaresult'.
#' @param annotation_file Character, the file path where the annotations will be written.
#' @param level Character, taxonomic level at which annotations are generated (default "Genus").
#'
#' @return None; this function writes annotations directly to the specified file.
#'
#' @details
#' The function starts by ensuring the `phyloseq` object is properly formatted with taxa as rows.
#' It calculates positive and negative counts for each taxon's representation in the results and
#' selects taxa based on a threshold. Each selected taxon's annotation includes visual attributes
#' like color and shape, tailored based on its statistical representation. These annotations are
#' designed to enhance visualizations, particularly phylogenetic trees, by providing clear indications
#' of the taxa's significance in the study.
#'
#' The function also checks for taxonomic rank format correctness and handles unknown values gracefully,
#' ensuring that the output file contains only well-formed and meaningful data. The colors and shapes
#' used in the annotations can be adjusted within the function to match specific visualization needs.
#'
#'
#' @importFrom phyloseq tax_table otu_table
#' @export
create_tax_selected <- function(ps, nburnin, niter, result, annotation_file, level="Genus"){
  if (is.null(level)){
    ps <- ps
  }else{
    ps <- tax_glom(ps,taxrank = level)
  }

  mytax <- as.data.frame(tax_table(ps))
  myotu <- as.data.frame(otu_table(ps))

  if (length(result$gammaresult)==nrow(myotu)){
    print("Your result and phyloseq subject are at the same level, continue to generate annotation file.")
  } else {
    stop("Please check your result and level, make sure they are at the same level.")
  }
  # ## calculate the negative counts for beta
  # betahat <- result$betahat[,nburnin:(nburnin+niter)]
  # rownames(betahat) <- rownames(myotu)
  # num_p <- apply(betahat,1,function(x) sum(x>0))
  # num_n <- apply(betahat, 1, function(x) sum(x<0))
  # df <- data.frame(row.names = rownames(myotu), positive=num_p, negative=num_n)
  # selected <- which(result$gammaresult>niter/2)
  #
  # otu_ann <- row.names(myotu[selected,])
  # otu_color <- df[rownames(df) %in% otu_ann,]
  # otu_color$color <- ifelse(otu_color$positive/niter>0.5, "red", "blue")
  otu_all <- row.names(myotu)
  ## create a data frame with color
  effect_size <- effect_size(result, ps, nburnin, niter, mode="mean",level=level)
  merged_df <- merge(effect_size, mytax,by=level,all.x = TRUE)
  merged_df$color <- ifelse(merged_df$effect_size > 0, "red", "blue")

  # Define the correctly formatted taxonomic rank names
  correct_format_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Actual ranks from the data frame
  actual_ranks <- colnames(mytax)

  # Identify ranks that are present but incorrectly formatted
  incorrectly_formatted_ranks <- actual_ranks[!actual_ranks %in% correct_format_ranks]

  # Generate an error message if there are any incorrectly formatted ranks
  if (length(incorrectly_formatted_ranks) > 0) {
    stop(paste("The following column names do not match the expected format:",
               paste(incorrectly_formatted_ranks, collapse=", "), ". Expected formats are:",
               paste(correct_format_ranks, collapse=", "), ". Please check your column names for correct capitalization."))
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
  ## write general annotation first
  writeLines("title\tTaxonomy tree with selection result", file_conn)
  writeLines("title_font_size\t25", file_conn)
  writeLines("start_rotation\t90", file_conn)
  # writeLines("clade_separation\t0.35", file_conn)
  writeLines("class_legend_font_size\t12", file_conn)
  writeLines("annotation_legend_font_size\t15", file_conn)
  writeLines("annotation_font_size\t10", file_conn)
  # writeLines("clade_marker_size\t0", file_conn)
  writeLines("branch_color_from_ancestor\t0", file_conn)
  # writeLines("annotation_background_separation\t-0.24", file_conn)
  # writeLines("annotation_background_offset\t-0.21", file_conn)
  # writeLines("annotation_background_width\t0.03", file_conn)
  writeLines("branch_thickness\t1.5", file_conn)
  writeLines("branch_bracket_depth\t0.5", file_conn)
  skip_count <-0
  for (otu in otu_all) {
    # species <- tax_tree$tax_paths[rownames(tax_tree)==otu]
    species <- as.character(mytax[otu,]$Species)
    class <- as.character(mytax[otu,]$Class)
    order <- as.character(mytax[otu,]$Order)
    family <- as.character(mytax[otu,]$Family)
    kingdom <- as.character(mytax[otu,]$Kingdom)
    level_tax <- as.character(mytax[otu,colnames(mytax)==level])
    # print(species)
    phylum <- as.character(mytax[otu,]$Phylum)
    genus <- as.character(mytax[otu,]$Genus)

    if (species == "unknown" | class == "unknown" | order == "unknown" | family == "unknown" |
        kingdom == "unknown" | phylum == "unknown" | genus == "unknown" | level_tax == "unknown") {
      skip_count <- skip_count+1
      next
    }

    color <- color_list$phylum_colors[color_list$unique_phyla==phylum]
    # print(otu)
    size <- ifelse(!is.na(level_abundance$size[level_abundance$expected_level == level_tax]),
                   level_abundance$size[level_abundance$expected_level == level_tax],
                   50)
    # phylum <- as.character(mytax[otu,]$Phylum)
    # family <- as.character(mytax[otu,]$Family)

    first_letter_genus <- toupper(substr(genus, 1, 1))
    first_letter_class <- toupper(substr(genus,1,1))
    first_letter_phylum <- toupper(substr(phylum,1,1))
    # phylum_class <- paste(phylum,class,sep = ".")
    # phylum_species <- paste(phylum, class,order,family,genus, species, sep = ".")
    # print(phylum_species)
    # phylum_genus <- paste(phylum, class, order, family,genus, sep = ".")
    # color <- phylum_colors[which(unique_phyla == phylum)]
    # writeLines(paste(species,"clade_marker_size", "50", sep = "\t"), file_conn)
    # writeLines(paste(kindom,"clade_marker_size", "0", sep = "\t"), file_conn)
    writeLines(paste(phylum,"clade_marker_size", "0", sep = "\t"), file_conn)
    writeLines(paste(order,"clade_marker_size", "0", sep = "\t"), file_conn)
    writeLines(paste(family,"clade_marker_size", "0", sep = "\t"), file_conn)
    # writeLines(paste(genus,"clade_marker_size", "0", sep = "\t"), file_conn)
    writeLines(paste(class,"clade_marker_size", "0", sep = "\t"), file_conn)
    writeLines(paste(genus, "ring_height", "1", size/10000, sep = "\t"), file_conn)

    writeLines(paste(genus,"ring_color", "1","purple", sep = "\t"), file_conn)
    writeLines(paste(genus,"clade_marker_color", "white", sep = "\t"), file_conn)

    writeLines(paste(phylum,"clade_marker_color","white", sep = "\t"), file_conn)
    writeLines(paste(phylum, "annotation_background_color","white", sep = "\t"),file_conn)
    # writeLines(paste(phylum,"annotation",paste(first_letter_phylum,":",phylum), sep = "\t"), file_conn)

    writeLines(paste(class, "clade_marker_color", "white", sep = "\t"), file_conn)

    # writeLines(paste(genus, "annotation", genus,sep = "\t"), file_conn)
    writeLines(paste(genus, "annotation_background_color","white", sep = "\t"),file_conn)
    writeLines(paste(genus, "annotation_rotation","90", sep = "\t"), file_conn)
    writeLines(paste(genus,"clade_marker_color", "white", sep = "\t"), file_conn)
  }
  print(paste("Total skips:", skip_count))
  ##write annotation for selected otu
  # for (otu in rownames(otu_color)) {
  #   species <- as.character(mytax[otu,]$Species)
  #   class <- as.character(mytax[otu,]$Class)
  #   order <- as.character(mytax[otu,]$Order)
  #   family <- as.character(mytax[otu,]$Family)
  #   kingdom <- as.character(mytax[otu,]$Kingdom)
  #   # print(species)
  #   phylum <- as.character(mytax[otu,]$Phylum)
  #   level_tax <- as.character(mytax[otu,colnames(mytax)==level])
  #   genus <- as.character(mytax[otu,]$Genus)
  #   if (species == "unknown" | class == "unknown" | order == "unknown" | family == "unknown" |
  #       kingdom == "unknown" | phylum == "unknown" | genus == "unknown" | level_tax == "unknown") {
  #     next
  #   }
  #   color1 <- otu_color$color[rownames(otu_color)==otu]
  #   # print(otu)
  #   size <- ifelse(!is.na(level_abundance$size[level_abundance$expected_level == level_tax]),
  #                  level_abundance$size[level_abundance$expected_level == level_tax],
  #                  50)
  #   # phylum <- as.character(mytax[otu,]$Phylum)
  #   # family <- as.character(mytax[otu,]$Family)
  #
  #   first_letter_genus <- toupper(substr(genus, 1, 1))
  #   first_letter_class <- toupper(substr(genus,1,1))
  #   first_letter_phylum <- toupper(substr(phylum,1,1))
  #   writeLines(paste(family,"clade_marker_shape","^", sep = "\t"), file_conn)
  #   writeLines(paste(family,"clade_marker_size","300", sep = "\t"), file_conn)
  #   writeLines(paste(family,"clade_marker_color", "darkgreen", sep = "\t"), file_conn)
  #   writeLines(paste(genus, "clade_marker_shape","^", sep="\t"), file_conn)
  #   writeLines(paste(genus, "annotation_background_color", "red", sep = "\t"), file_conn)
  #   writeLines(paste(class, "annotation_background_color", "red", sep = "\t"), file_conn)
  #   writeLines(paste(genus, "clade_marker_size", "300", sep = "\t"), file_conn)
  #   writeLines(paste(genus,"clade_marker_color", "darkgreen", sep = "\t"), file_conn)
  #   # writeLines(paste(genus, "annotation", paste("Genus_",first_letter_genus,":",genus),sep = "\t"), file_conn)
  #   writeLines(paste(class, "clade_marker_shape", "^", sep = "\t"), file_conn)
  #   # writeLines(paste(class, "annotation_background_color", "red", sep = "\t"), file_conn)
  #   writeLines(paste(class, "clade_marker_size", "300", sep = "\t"), file_conn)
  #   writeLines(paste(class,"clade_marker_color", "darkgreen", sep = "\t"), file_conn)
  #   writeLines(paste(phylum, "clade_marker_shape","^", sep = "\t"), file_conn)
  #   # writeLines(paste(phylum, "annotation_background_color", "red", sep = "\t"), file_conn)
  #   writeLines(paste(phylum, "clade_marker_size", "300", sep = "\t"), file_conn)
  #   writeLines(paste(phylum,"clade_marker_color", "darkgreen", sep = "\t"), file_conn)
  #   writeLines(paste(phylum, "annotation_background_color", "red", sep = "\t"), file_conn)
  #   writeLines(paste(phylum,"annotation",paste("Phylum_",first_letter_phylum,":",phylum), sep = "\t"), file_conn)
  #
  # }
  for (i in 1:nrow(merged_df)) {
    # Extract taxonomic information for the current OTU
    otu <- merged_df$level[i]
    species <- merged_df$Species[i]
    class <- merged_df$Class[i]
    order <- merged_df$Order[i]
    family <- merged_df$Family[i]
    kingdom <- merged_df$Kingdom[i]
    phylum <- merged_df$Phylum[i]
    genus <- merged_df$Genus[i]
    effect_size <- merged_df$effect_size[i]

    # Skip if any taxonomic level is "unknown"
    if ("unknown" %in% c(species, class, order, family, kingdom, phylum, genus)) {
      next
    }

    # Determine the color based on effect size
    color <- ifelse(effect_size > 0, "red", "blue")

    ## write annotation
    first_letter_phylum <- toupper(substr(phylum,1,1))
    writeLines(paste(otu, "clade_marker_shape","^", sep="\t"), file_conn)
    writeLines(paste(otu, "annotation_background_color", "yellow", sep = "\t"), file_conn)
    writeLines(paste(otu, "clade_marker_size", "300", sep = "\t"), file_conn)
    writeLines(paste(otu,"clade_marker_color", color, sep = "\t"), file_conn)
    # writeLines(paste(genus, "annotation", paste("Genus_",first_letter_genus,":",genus),sep = "\t"), file_conn)
    writeLines(paste(phylum, "clade_marker_shape","^", sep = "\t"), file_conn)
    # writeLines(paste(phylum, "annotation_background_color", "red", sep = "\t"), file_conn)
    writeLines(paste(phylum, "clade_marker_size", "300", sep = "\t"), file_conn)
    writeLines(paste(phylum,"clade_marker_color", "darkgreen", sep = "\t"), file_conn)
    writeLines(paste(phylum, "annotation_background_color", "red", sep = "\t"), file_conn)
    writeLines(paste(phylum,"annotation",paste("Phylum_",first_letter_phylum,":",phylum), sep = "\t"), file_conn)
  }
  close(file_conn)
}
