#' Bayesian-tree based Model Fitting for Microbiome Data
#'
#' This function processes microbiome data, aggregates it at specified taxonomic levels,
#' and fits a Bayesian-tree based model. It checks the integrity of the
#' operational taxonomic unit (OTU) table, scales the data, prepares the phylogenetic matrix,
#' and runs a Gibbs sampling algorithm to estimate model parameters.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#' @param nburnin Integer, number of burn-in iterations for the Gibbs sampler.
#' @param niter Integer, number of sampling iterations after burn-in.
#' @param a Numeric vector or scalar, hyperparameters for the model (choose it carefully!).
#' @param level Optional character, taxonomic level at which to aggregate data (e.g., "Phylum").
#'        If NULL, uses the original taxonomy without aggregation.
#' @param response Character, name of the response variable in the sample data.
#'
#' @return Returns a list containing the results from the Gibbs sampler, including parameter estimates
#'         and diagnostics from the sampling process.
#'
#' @details
#' The function first validates the OTU table structure, ensuring that taxa are properly represented.
#' It then optionally aggregates the OTU data to a specified taxonomic level. Response variables and
#' predictors are extracted from the `phyloseq` object. Data is then log-transformed and standardized.
#' A phylogenetic variance-covariance matrix is computed and adjusted. The function finally invokes the
#' `gibbsgamma` function for fitting the Bayesian model.
#'
#'
#' @importFrom phyloseq sample_data otu_table taxa_names phy_tree tax_glom transform_sample_counts otu_table<-
#' @import Matrix
#' @importFrom ape vcv
#' @export
result <-
  function(ps,nburnin,niter,a,level = NULL,response){
    ## check whether otu table is correct
    if (!all(taxa_names(ps) %in% rownames(otu_table(ps)))) {
      stop("Error: Taxa are not represented as rows in the data.")
    }
    myotu <- as.data.frame(otu_table(ps))
    myotu[myotu==0]=0.5
    myotu <- myotu/rowSums(myotu)
    otu_table(ps) <- otu_table(myotu,taxa_are_rows = TRUE)
    if (is.null(level)){
      ps1 <- ps
    }else{
      ps1 <- tax_glom(ps,taxrank = level)
    }
      seed <- 101

      ## Obtain response variable and predictors
      myotu =as.data.frame(otu_table(ps1))
      mysam = as.data.frame(sample_data(ps1))
      Z = t(myotu)
      Y = mysam[[response]]

      Z = as.matrix(Z)
      stand <- list()
      Z <- scale(log(Z))
      stand$mux <- attr(Z, "scaled:center")
      stand$Sx <- attr(Z, "scaled:scale")
      X <- Z
      z <- as.matrix(Y)
      z <- scale(z)
      stand$muy <- attr(z, "scaled:center")
      stand$Sy <- attr(z, "scaled:scale")
      Y <- z

      n <- nrow(X)
      p <- ncol(X)

      c = 100
      N <- rbind(diag(p), rep(c, p)) %*% diag(1/stand$Sx)

      nop <- floor(min(n,p)/2)

      ## prepare phylogentic matrix
      Q1 <- vcv(phy_tree(ps1))
      # Q1 <- read.csv("sortcorr.csv", header = F, sep = "\t")

      Q1[Q1 > 0.9] <- 0.9
      Q1 <- Q1 + 1 * diag(p)
      Q1 <- as.matrix(Q1)
      Q1 <- solve(Q1)
      Q1 <- Q1 - diag(diag(Q1))
      Q1[Q1 > 1] = 1
      Q1[Q1 < -1] = -1
      # sqrt(sum(Q1>0))
      Q <- as(Q1, "sparseMatrix")

      tau=1
      nu=0
      omega=0

      ss = 50
      al=-15
      au=0
      # a0 = seq(from = al, to = au, by = (au-al) / ss)
      a <- rep(a,p)

      result <- gibbsgamma(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, TRUE, stand, TRUE)
      return(result)
  }
