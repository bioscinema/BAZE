#' Aggregate and Visualize Selection Probabilities
#'
#' This function aggregates results from a Gibbs sampling process and visualizes selection
#' probabilities of taxa at specified taxonomic levels, comparing positive and negative selections.
#' It calculates the counts of negative and positive selections for each taxon, merges these with
#' taxonomic data, and produces a bar plot grouped by the specified taxonomic level.
#'
#' @param nburnin Integer, number of burn-in iterations to skip in the analysis.
#' @param niter Integer, number of iterations considered for aggregating results.
#' @param result List, containing results from a Gibbs sampling algorithm, specifically `betahat`
#'        for estimated coefficients and `gammaresult` for variable inclusion results.
#' @param taxa_table Data frame, containing taxonomic information for each taxon. Must include
#'        a column corresponding to `level` that contains taxonomic classifications.
#' @param level Character, specifying the taxonomic level to aggregate and visualize data.
#'        Defaults to "Phylum".
#' @param label_size Numeric, size of the label text in the plot, default is 10.
#'
#' @return A ggplot object displaying the aggregated selection probabilities as a bar plot,
#'         grouped by the specified taxonomic level. Each group is annotated with the number of
#'         taxa included in that group.
#'
#' @examples
#' \dontrun{
#'   result <- list(betahat = matrix(rnorm(2000), 100, 20),
#'                  gammaresult = sample(0:1, 100, replace = TRUE))
#'   taxa_table <- data.frame(Phylum = rep(c("Firmicutes", "Bacteroidetes"), each = 50))
#'   plot <- aggregate_plot(nburnin = 500, niter = 1500, result = result,
#'                          taxa_table = taxa_table, level = "Phylum")
#'   print(plot)
#' }
#' @import ggplot2
#' @import dplyr
#' @export
aggregate_plot <- function(nburnin, niter , result, taxa_table, level = "Phylum", label_size=10){
  ## calculate the negative counts for beta
  betahat <- result$betahat[,nburnin:(nburnin+niter)]
  betahat_n <- which(betahat<0, arr.ind = TRUE)
  counts_n <- as.data.frame(table(betahat_n[,1]))
  colnames(counts_n) <- c("otu_id", "Freq")
  ## calculate the total counts for beta
  gamma <- result$gammaresult
  gamma <- as.data.frame(which(gamma > 0))
  colnames(gamma) <- c("otu_id")
  gamma$betahat <- result$gammaresult[which(result$gammaresult > 0)]
  ## merge two in one table
  betahat_c <- merge(gamma, counts_n, by = "otu_id", all.x = TRUE)
  betahat_c[is.na(betahat_c)] <- 0
  colnames(betahat_c)[3] <- c("neg_counts")
  # betahat_c <- betahat_c %>%
  #             rename(
  #               neg_counts = Freq
  #             )
  row.names(betahat_c) <- betahat_c$otu_id
  row.names(taxa_table) <- 1:nrow(taxa_table)
  new_tax <- merge(betahat_c, taxa_table, by = "row.names", all.x = TRUE)
  ## create data frame for bar plot
  if (! level %in% colnames(new_tax)) {
    stop(paste("error:", level, "not found in your data frame"))
  }
  else {
    new_tax <- new_tax[order(new_tax[[level]]),]
    group <- rep(new_tax[[level]], each = 2)
    tem <- data.frame(group)
    barc <- tem %>%
      group_by(group) %>%
      mutate(counter = rep(1:ceiling(n()/2), each = 2, length.out = n()))
    bar <- barc$counter
    value <- rep(0, length(group))
    for (i in 1:length(new_tax[[level]])) {
      p = 2*i - 1
      n = 2*i
      value[p] = new_tax$betahat[i] / (niter + 1) * (1 - new_tax$neg_count[i] / new_tax$betahat[i])
      value[n] =new_tax$betahat[i] / (niter + 1) * (new_tax$neg_count[i] / new_tax$betahat[i])
    }
    segment <- factor(rep(c("positive", "negative"), times = length(new_tax$otu_id)), levels = c("positive", "negative"))
    df <- data.frame(group, bar, value, segment)
    group_counts <- table(df$group)
    df_summary <- df %>%
      group_by(group) %>%
      summarise(counts = n() / 2)
    df$group <- factor(df$group, levels = df_summary$group, labels = paste(df_summary$group, "(counts=", df_summary$counts, ")", sep = ""))
  }
  p <- ggplot(df, aes(x = factor(bar), y = value, fill = segment, group = group)) +
    geom_bar(stat = "identity", position = "stack", width = 0.4) +
    scale_fill_manual(values = c("lightblue", "lightcoral")) +
    theme_minimal() +
    labs(title = "selection prob plot", y = "selection prob", x = "Group") +
    facet_wrap(~ group, scales = "free_x", strip.position = "bottom") +
    scale_x_discrete("Groups",
                     breaks = unique(interaction(df$group)),
                     labels = rep(paste(df_summary$group, "(counts=", df_summary$counts, ")", sep = ""), each = 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length = unit(0.25, "cm"),
          strip.text.x = element_text(size=label_size))
  print(p)
}
