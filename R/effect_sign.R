#' Generate Effect Significance for Taxonomic Levels
#'
#' This function computes the effect size and significance (positive or negative) for taxa at a specified taxonomic level based on the provided phyloseq object and model results.
#'
#' @param result A list containing the `betahat` matrix and `gammaresult` vector from the model results.
#' @param ps A phyloseq object containing the microbiome data.
#' @param nburnin An integer specifying the number of burn-in iterations.
#' @param niter An integer specifying the total number of iterations.
#' @param mode A character string specifying the mode of calculation (currently unused, reserved for future use). Default is "mean".
#' @param level A character string specifying the taxonomic level at which to aggregate the data (e.g., "Genus"). If NULL, no aggregation is performed.
#' @return A data frame containing taxa names and their corresponding effect significance ("positive" or "negative").
#' @details
#' This function first aggregates the taxa in the phyloseq object to the specified taxonomic level (if provided). It then ensures the dimensions of the taxonomy table align with the dimensions of the `betahat` matrix from the model results. The function combines the taxa names with the selection results (`gammaresult`), filters taxa based on a threshold, calculates the effect size for each taxon, and determines the effect's significance (positive or negative). Finally, it returns a data frame with taxa names and their effect significance, excluding taxa labeled as "unknown" or "uncultured".
#' @examples
#' # Assuming 'result' is a list with 'betahat' and 'gammaresult', 'ps' is a phyloseq object
#' # effect_sign(result, ps, nburnin=100, niter=1000, level="Genus")
#' @importFrom stats na.omit
#' @export
effect_sign <- function(result, ps, nburnin, niter, mode = "mean", level = "Genus") {
  # Aggregate taxa to the specified taxonomic level if provided
  if (is.null(level)) {
    ps1 <- ps
  } else {
    ps1 <- tax_glom(ps, taxrank = level)
  }
  
  # Extract taxonomy table and ensure alignment with result dimensions
  mytax <- as.data.frame(tax_table(ps1))
  if (nrow(result$betahat) != nrow(mytax)) {
    stop("Please re-enter level to keep your result and phyloseq subject aligned")
  }
  
  # Extract betahat and gammaresult from result
  betahat <- result$betahat
  gammaresult <- result$gammaresult
  
  # Combine taxa name with selection result
  gamma <- data.frame(taxa = mytax[[level]], result = gammaresult)
  combined_result <- cbind(gamma, betahat)
  
  # Filter results based on gamma threshold
  betahat_s <- combined_result[combined_result$result > niter / 2, ]
  
  # Calculate effect size
  effect <- rowMeans(betahat_s[, (nburnin + 3):ncol(betahat_s)])
  betahat_s$sign <- ifelse(effect > 0, "positive", "negative")
  
  # Generate effect size data frame and filter out unwanted taxa
  effect_sign <- data.frame(taxa = betahat_s$taxa, effect_sign = betahat_s$sign)
  names(effect_sign) <- c(level, "effect_sign")
  effect_sign <- effect_sign[!effect_sign[[level]] %in% c("unknown", "uncultured"), ]
  
  # Return the result
  return(effect_sign)
}