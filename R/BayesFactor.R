#' Calculate Bayes Factor for Model Comparison
#'
#' This function computes the Bayes Factor for comparing models in the context of Gibbs sampling
#' for Bayesian variable selection. It supports handling updates to the model based on the inclusion
#' or exclusion of a variable (flagged by the `flag` parameter). The function also updates
#' and returns the auxiliary information used in the Gibbs sampler.
#'
#' @param Y Numeric vector or matrix, the response variables.
#' @param Xri Numeric matrix, the design matrix without the variable of interest.
#' @param Xi Numeric vector or matrix, the design matrix for the variable of interest.
#' @param XIi Numeric matrix, the combined design matrix with the variable of interest.
#' @param Tri Numeric matrix, transformation matrix associated with Xri.
#' @param Ti Numeric vector or matrix, transformation matrix associated with Xi.
#' @param TIi Numeric matrix, combined transformation matrix.
#' @param invAri Numeric matrix, the inverse of the covariance matrix of the regression coefficients without the variable of interest.
#' @param n Integer, sample size or number of observations.
#' @param tau Numeric, hyperparameter scaling the precision matrix.
#' @param nu Numeric, degrees of freedom for the prior distribution of the error variance.
#' @param omega Numeric, scale parameter for the prior distribution of the error variance.
#' @param flag Boolean, indicates whether to include or exclude the variable of interest in the model.
#' @param keep List, contains previously calculated values required for updates (residual sum of squares and determinants).
#'
#' @return A list containing:
#'   - `F1` or `F`: The calculated Bayes Factor, depending on the flag.
#'   - `keep`: Updated list of auxiliary information.
#'   - `invAi`: Optionally, the updated inverse covariance matrix (only if `flag` is FALSE).
#'
#' @details
#' The function performs a rank-1 update to the Cholesky decomposition when `flag` is FALSE,
#' indicating a variable exclusion scenario. When `flag` is TRUE, it simply recalculates based on
#' existing values. This functionality is crucial for dynamic model updating in Bayesian variable
#' selection techniques.
#'
#' @export
BayesFactor <- function(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep) {
  keep <- keep
  # lri1 <- Rfast::cholesky(invAri, parallel = TRUE)
  # Assuming invAri is your matrix
  lri1 <- chol(invAri)
  # print(max(eigen(lri1)$values))
  Lri <- t(lri1)

  sqrtdetinvAri <- sum(log(diag(Lri)))
  resAri <- t(Y) %*% Y - t(Y) %*% Xri %*% invAri %*% t(Xri) %*% Y

  if (flag) {
    # copy invAi from previous step
    # keep1 <- keep
    sqrtdetinvAi <- keep$sqrtdetinvAi
    resAi <- keep$resAi

    keep$resAri <- resAri
    keep$sqrtdetinvAri <- sqrtdetinvAri

    logratiodetT <- sum(log(diag(chol(t(Tri) %*% Tri)))) - sum(log(diag(chol(t(TIi) %*% TIi))))

    # invAi <- invAi
    F1 <- -log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT + (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi))
    F1 <- exp(F1)
    return(list(F1 = F1, keep = keep))
  } else {
    # update invAi from invAri
    # keep <- keep
    Srii <- t(Xri) %*% Xi + tau^(-2) * (t(Tri) %*% Ti)
    sii <- t(Xi) %*% Xi + tau^(-2) * (t(Ti) %*% Ti)

    # calculate inverse Ai^(-1) given Ari^(-1)
    v1 <- sqrt(1 / (sii * (1 - t(Srii) %*% invAri %*% Srii / sii)))
    v <- v1[1,1] * invAri %*% Srii
    #v <- sqrt(1 / (sii * (1 - t(Srii) %*% invAri %*% Srii / sii))) %*% invAri %*% Srii

    A11 <- invAri + v %*% t(v)
    A12 <- -A11 %*% (Srii / sii[1,1])
    A21 <- -(1 / sii) %*% t(Srii) %*% A11
    A22 <- 1 / sii[1,1] + (1 / sii[1,1]) %*% t(Srii) %*% A11 %*% Srii / sii[1,1]
    invAi <- rbind(cbind(A11, A12), cbind(A21, A22))
    resAi <- t(Y) %*% Y - t(Y) %*% XIi %*% invAi %*% t(XIi) %*% Y

    # calculate determinant of Ai^(-1) and Ari^(-1)
    # Rank 1 update to Cholesky factorization
    tilLri <- t(chol(Lri %*% t(Lri)+v %*% t(v)))
    # print(tilLri)
    Lrii <- forwardsolve(tilLri, A12)
    # Lrii <- chol2inv(tilLri, A12)
    lii <- sqrt(A22 - t(Lrii) %*% Lrii)
    Li <- rbind(cbind(tilLri, matrix(0, nrow(tilLri), 1)), cbind(t(Lrii), lii))
    # print(Li)
    sqrtdetinvAi <- sum(log(diag(Li)))

    # keep quantities for future steps
    keep$resAri <- resAri
    keep$sqrtdetinvAri <- sqrtdetinvAri

    keep$resAi <- resAi
    keep$sqrtdetinvAi <- sqrtdetinvAi
    logratiodetT <- -sum(log(diag(chol(t(Tri) %*% Tri)))) + sum(log(diag(chol(t(TIi) %*% TIi))))

    F <- -log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT + (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi))
    F <- exp(F)
    return(list(F = F, keep = keep, invAi = invAi))
  }

  # calculate the Bayes factor: bf = P(Y|gammai=0,gamma_-i) / P(Y|gammai=1,gamma_-i)
  # F <- -log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT + (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi))
  # F <- exp(F)

  # return(list(F = F, keep = keep, invAi = invAi))
}
