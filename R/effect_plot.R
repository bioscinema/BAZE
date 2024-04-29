#' Plot Effect Sizes from Bayesian Model Results
#'
#' This function calculates and plots the effect sizes derived from Bayesian variable
#' selection model results. It processes the results to compute effect sizes based on
#' specified summary statistics (mean or median) after aggregating data at a given taxonomic level.
#' The effect sizes are then visualized in a bar plot with bars colored based on their
#' sign (positive or negative).
#'
#' @param result A list from Bayesian variable selection containing at least the betahat
#'        coefficients and the gamma result selection frequencies.
#' @param ps A `phyloseq` object containing microbiome data which will be aggregated
#'        at the specified taxonomic level.
#' @param nburnin The number of burn-in iterations in the Bayesian model; these initial iterations
#'        will be excluded from the effect size calculation.
#' @param niter The total number of iterations after burn-in, used to determine the threshold
#'        for significant selection in the visualization.
#' @param mode A character string specifying the method to summarize the effect sizes across iterations.
#'        Options are "mean" or "median". Defaults to "mean".
#' @param level The taxonomic level at which to aggregate data in the `phyloseq` object.
#'        Defaults to "Genus".
#'
#' @return A ggplot object representing the effect size plot, with taxa on the y-axis and
#'         effect sizes on the x-axis. Bars are filled according to the sign of the effect size.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual labs theme_minimal theme
#' @importFrom stats reorder
#' @export
effect_plot <- function(result, ps, nburnin, niter, mode="mean",level="Genus"){
  ##calculate effect size
  effect_size <- effect_size(result, ps, nburnin, niter, mode=mode,level=level)

  ##add a sign column
  effect_size$Sign <- ifelse(effect_size$effect_size > 0, "Positive", "Negative")


  ##generate effect size plot
  p <- ggplot(effect_size, aes(x = reorder(level, effect_size), y = effect_size, fill=Sign)) +
    geom_col(color="black", width = 0.7) +  # Using geom_col which is geom_bar(stat = "identity")
    coord_flip() +  # Flip coordinates for horizontal bars
    # scale_fill_manual(values = c("Negative" = "red", "Positive" = "blue"),
    #                   name = "Effect Sign",
    #                   labels = c("Negative", "Positive")) +
    labs(x = "Taxa", y = "Effect Size", title = "Effect Size Plot") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ##return plot
  return(p)
}
