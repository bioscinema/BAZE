#' Calculate Effect Sizes for Bayesian Variable Selection Model Results
#'
#' This function calculates effect sizes for variables in a Bayesian variable selection model
#' based on the `betahat` coefficients and selection frequency across iterations. The function
#' allows users to choose between mean or median to summarize the effect sizes and requires
#' that the number of taxa in the `phyloseq` object matches the number of `betahat` coefficients.
#'
#' @param result A list containing the results from a Bayesian variable selection model,
#'   expected to include elements `betahat` and `gammaresult` where `betahat` are the
#'   coefficient estimates across iterations and `gammaresult` indicates selection frequency.
#' @param ps A `phyloseq` object containing microbiome data, which will be aggregated at the specified
#'   taxonomic level before analysis.
#' @param nburnin Integer specifying the number of burn-in iterations in the model, used to
#'   exclude these iterations from effect size calculations.
#' @param niter Integer specifying the total number of iterations after burn-in, used to determine
#'   the threshold for considering a variable significantly selected.
#' @param mode Character string specifying how to summarize the effect sizes across iterations.
#'   Options are "mean" or "median". Defaults to "mean".
#' @param level The taxonomic level at which to aggregate data in the `phyloseq` object.
#'   Defaults to "Genus".
#'
#' @return Returns a data frame with taxa and their corresponding effect sizes, ordered
#'   by effect size in decreasing order.
#'
#' @importFrom phyloseq tax_glom tax_table
#' @importFrom stats median
#' @export
effect_size <- function(result, ps, nburnin,niter, mode="mean",level="Genus"){
  if (is.null(level)){
    ps1 <- ps
  }else{
    ps1 <- tax_glom(ps,taxrank = level)
  }
  mytax <- as.data.frame(tax_table(ps1))
  if (nrow(result$betahat) != nrow(mytax)) {
    stop("please re-enter level to keep your result and phyloseq subject align")
  }

  ## extract betahat and frequencies from result
  betahat <- result$betahat
  gammaresult <- result$gammaresult

  ##combine taxa name with selection result
  gamma <- data.frame(taxa = mytax[[level]], result = gammaresult)
  result <- cbind(gamma,betahat)
  betahat_s <- result[which(result$result>niter/2),]

  ##calculate effect size
  if (mode == "mean") {
    betahat_s$effect = rowMeans(betahat_s[,(nburnin+3):ncol(betahat_s)])
  } else if(mode=="median"){
    betahat_s$effect <- apply(betahat_s[,(nburnin+3):ncol(betahat_s)], 1, median)
  }

  ##generate effect size data frame and order it
  effect_size <- data.frame(taxa=betahat_s$taxa,effect_size=betahat_s$effect)
  names(effect_size) <- c(level, "effect_size")
  effect_size <- effect_size[-effect_size[[level]] %in% c("unknown","uncultured"),]
  effect_size <- effect_size[order(effect_size$effect_size,decreasing = TRUE),]

  ##return the result
  return(effect_size)
}
