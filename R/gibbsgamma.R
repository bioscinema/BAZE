#' Gibbs Sampling for Bayesian Variable Selection
#'
#' This function implements a Gibbs sampling algorithm for variable selection in a Bayesian
#' framework, using a gamma-shrinkage prior. It splits the data into training and testing sets,
#' performs the sampling, and optionally predicts responses based on selected models.
#'
#' @param nburnin Integer, number of burn-in iterations for the Gibbs sampler.
#' @param niter Integer, number of iterations after burn-in.
#' @param p Integer, number of predictors.
#' @param nop Integer, initial number of predictors to include.
#' @param Y Numeric vector or matrix, response variable(s).
#' @param X Numeric matrix, predictor variables.
#' @param N Numeric matrix, normalization matrix for predictors.
#' @param a Numeric vector, hyperparameters for the gamma prior (usually negative number).
#' @param Q Numeric matrix, phylogenetic matrix.
#' @param n Integer, number of observations in Y and rows in X.
#' @param tau Numeric, hyperparameter scaling the precision matrix N.
#' @param nu Numeric, degrees of freedom for the inverse gamma distribution of sigma^2.
#' @param omega Numeric, scale parameter for the inverse gamma distribution of sigma^2.
#' @param seed Integer, seed for random number generation to ensure reproducibility.
#' @param predict Logical, whether to perform prediction on the test set.
#' @param stand List, contains standardized parameters including 'mux' and 'Sx' for X, and 'muy' and 'Sy' for Y.
#' @param display Logical, whether to display progress of the Gibbs sampler.
#' @param split_rate Numeric, proportion of the data to be used as the training set (default is 0.7).
#'
#' @importFrom stats rbinom
#' @return A list containing:
#'   - `gamma`: Matrix of inclusion indicators for variables across all iterations.
#'   - `betahat`: Matrix of estimated coefficients.
#'   - `gammaresult`: Sum of gamma indicators across post-burn-in iterations.
#'   - `MSE`: Mean squared error of predictions on the testing set.
#'   - `nselect`: Number of selected variables in each iteration.
#'   - `Yhat`: Predicted values for the testing set.
#' @export
gibbsgamma <- function(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, predict, stand, display, split_rate=0.7) {

  # Setting the seed
  set.seed(seed)
  # nop=12
  # Determine indices for splitting the data (e.g., 70% training, 30% testing)
  train_indices <- sample(1:nrow(X), size = floor(split_rate * n), replace = FALSE)
  test_indices <- setdiff(1:nrow(X), train_indices)

  # Split the data into training and testing sets
  X_train <- X[train_indices, ]
  Y_train <- Y[train_indices,]

  X_test <- X[test_indices, ]
  Y_test <- Y[test_indices,]


  # Initialize index, nop is the number of initial 1's
  # Sample with replacement, p more than nop

  index <- sample(1:p, nop, replace = FALSE)


  index <- sort(index)
  nop <- length(index)

  index
  # Initialize gamma
  gamma <- matrix(0, nrow = nburnin + niter + 1, ncol = p)


  gamma[1, index] <- 1


  if (is.null(stand)) {
    invSx <- diag(p)
    Sy <- 1
    Yobs_train <- Y_train
    Yobs_test <- Y_test
    # Xobs <- X
  } else {
    invSx <- diag(1 / stand$Sx)
    Sy <- stand$Sy
    Yobs_train <- Y_train * stand$Sy + stand$muy
    Yobs_test <- Y_test * stand$Sy + stand$muy
    # Xobs <- X * diag(stand$Sx) + stand$mux
  }

  # Generate a sequence of index for update
  # set.seed(seed)

  proposeindx <- sample(1:p, nburnin + niter + 1, replace = TRUE)

  # Initialize Ari and keep
  Xri <- X_train[, index]
  Xri_test <- X_test[ ,index]
  Tri <- N[, index]
  Ari <- t(Xri) %*% Xri + tau^(-2) * (t(Tri) %*% Tri)
  Ari_test <- t(Xri_test) %*% Xri_test + tau^(-2) * (t(Tri) %*% Tri)
  invAri<-chol2inv(chol(Ari))
  # invAri <- solve(Ari)
  # print(dim(invAri))
  invAri_test <- chol2inv(chol(Ari_test))
  # print(dim(invAri_test))
 # invAri <- solve(Ari,diag(nop))
  Lri <- t(chol(invAri))

  keep <- list()
  keep$resAi <- t(Y_train) %*% Y_train - t(Y_train) %*% Xri %*% invAri %*% t(Xri) %*% Y_train
  keep$sqrtdetinvAi <- sum(log(diag(Lri)))

  if (predict) {
    # initialize beta
    betahat <- matrix(0, nrow = p, ncol = nburnin + niter + 1)
    tembeta <- invAri %*% t(Xri) %*% Y_train
    # print(dim(tembeta))
    # print(length(index))
    # print(dim(invSx[index,index]))
    betahat[index, 1] <- Sy * invSx[index, index] %*% tembeta
    Yhat <- matrix(0, nrow = length(test_indices), ncol = nburnin + niter + 1)
    # print(dim(stand$muy + Sy * X_test[, index] %*% tembeta))
    # print(length(Yhat[,1]))
    Yhat[, 1] <- stand$muy + Sy * X_test[, index] %*% tembeta
    MSE <- numeric(nburnin + niter + 1)
    MSE[1] <- 1 / n * t(Yobs_test - Yhat[, 1]) %*% (Yobs_test - Yhat[, 1])
    nselect <- numeric(nburnin + niter + 1)
    nselect[1] <- nop
  }


  if (display) {
    k <- 1
    cat("Gibbs Sampling is starting and will print every 500 iterations.\n")
  }

  for(i in 1:(nburnin+niter)){
  # for(i in 1:5){
    #print(i)
    gamma[i+1,] <- gamma[i,]
    betahat[,i+1] <- rep(0, p)
    flag <- any(index == proposeindx[i])
    if (flag) {
      ## proposed variable is inside existing index
      indxtemp <- index[index != proposeindx[i]]
      # print(indxtemp)
      Xri <- X_train[, indxtemp]
      Xi <- X_train[, proposeindx[i]]
      XIi <- cbind(Xri, Xi)
      Tri <- N[, indxtemp]
      Ti <- N[, proposeindx[i]]
      TIi <- cbind(Tri, Ti)
      invAritemp <- invAri
      idx <- list()
      tn <- length(index)
      #print(tn)
      seq <- 1:tn
      if (proposeindx[i] == max(index)) {
        idx[[1]] <- seq
        idx[[2]] <- seq
      } else if (proposeindx[i] == min(index)) {
        idx[[1]] <- c(2:tn, 1)
        idx[[2]] <- c(2:tn, 1)
      } else {
        ti <- seq[index == proposeindx[i]]
        idx[[1]] <- c(1:(ti - 1), (ti+1):tn, ti)
        idx[[2]] <- c(1:(ti - 1), (ti+1):tn, ti)
      }

      invAitemp <- invAri[idx[[1]], idx[[2]]]
      # swap the ti variable to the last one of the list
      invAri1 <- invAitemp[1:(tn-1), 1:(tn-1)]

      # print(dim(invAri1)[1])
      invAri2 <- invAitemp[1:(tn-1), tn, drop = FALSE]
      invAri3 <- invAitemp[tn, tn]
      invAri <- invAri1 - invAri2 %*% t(invAri2) / invAri3
      # invAri_1 <- invAritemp[1:(tn-1), 1:(tn-1)] - invAritemp[1:(tn-1), tn] %*% invAritemp[tn, 1:(tn-1)] / invAritemp[tn, tn]
      result1 <- BayesFactor(Y_train, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep)
      F1 <- result1$F1
      keep <- result1$keep
      pgammai1 <- exp(a[proposeindx[i]] + Q[proposeindx[i], indxtemp] %*% gamma[i, indxtemp])
      pcond <- 1 / (1 + 1 / (F1 * pgammai1))
      # print(pcond)
      newgamma <- rbinom(1, 1, pcond)
      gamma[i+1, proposeindx[i]] <- newgamma
      if ((newgamma == 0) && (tn > 1)) {
        index <- indxtemp
        if(length(index)==1){
          warning("Only 1 variable is left. Please increase the penalty value a.")
          break;
        }
        keep$resAi <- keep$resAri
        keep$sqrtdetinvAi <- keep$sqrtdetinvAri
        if (predict) {
          tembeta <- invAri %*% t(Xri) %*% Y_train
          # print(dim(tembeta))
          # print(length(index))
          # print(dim(invSx[index+1, index+1]))
          betahat[index, i + 1] <- Sy * invSx[index, index] %*% tembeta
          # print(dim(stand$muy + Sy * X_test[, index] %*% tembeta))
          # print(dim(Yhat[,i+1]))
          Yhat[, i + 1] <- stand$muy + Sy * X_test[, index] %*% tembeta
          MSE[i + 1] <- 1 / (n * t(Yobs_test - Yhat[, i + 1]) %*% (Yobs_test - Yhat[, i + 1])*(1-split_rate))
          nselect[i + 1] <- tn - 1
        }

      }else if ((newgamma != 0) || (tn == 1)) {
        invAri <- invAritemp
        if (predict) {
          betahat[, i + 1] <- betahat[, i]
          Yhat[, i + 1] <- Yhat[, i]
          MSE[i + 1] <- MSE[i]
          nselect[i + 1] <- nselect[i]
        }
      }


    }else {
      indxtemp <- index
      Xri <- X_train[, indxtemp]
      Xi <- X_train[, proposeindx[i]]
      XIi <- cbind(Xri, Xi)
      Tri <- N[, indxtemp]
      Ti <- N[, proposeindx[i]]
      TIi <- cbind(Tri, Ti)
      BayesResult <- BayesFactor(Y_train, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep)
      F <- BayesResult$F
      #print(F)
      keep <- BayesResult$keep
      invAi <- BayesResult$invAi
      pgammai1 <- exp(a[proposeindx[i]] + Q[proposeindx[i], indxtemp] %*% gamma[i, indxtemp])
      pcond <- 1 / (1 + 1/(F * pgammai1))
      #print(pcond)
      newgamma <- rbinom(1, 1, pcond)
      gamma[i + 1, proposeindx[i]] <- newgamma
      if (newgamma == 1) {
        #include the proposed variable
        index <- sort(c(indxtemp, proposeindx[i]))
        invAri <- invAi
        idx <- list()
        tn <- length(index)
        seq <- 1:tn
        if (proposeindx[i] > max(indxtemp)) {
          idx[[1]] <- seq
          idx[[2]] <- seq
        } else if (proposeindx[i] < min(indxtemp)) {
          idx[[1]] <- c(tn, 1:(tn - 1))
          idx[[2]] <- c(tn, 1:(tn - 1))
        } else {
          ti <- seq[index == proposeindx[i]]
          idx[[1]] <- c(1:(ti - 1), tn, ti:(tn - 1))
          idx[[2]] <- c(1:(ti - 1), tn, ti:(tn - 1))
        }
        invAri <- invAri[idx[[1]], idx[[2]]]
        if (predict) {
          tembeta <- invAri %*% t(X_train[, index]) %*% Y_train
          betahat[index, i + 1] <- Sy * invSx[index, index] %*% tembeta
          Yhat[, i + 1] <- stand$muy + Sy * X_test[, index] %*% tembeta
          MSE[i + 1] <- 1 / (n * t(Yobs_test - Yhat[, i + 1]) %*% (Yobs_test - Yhat[, i + 1])*(1-split_rate))
          nselect[i + 1] <- tn
        }
      }else {
        keep$resAi <- keep$resAri
        keep$sqrtdetinvAi <- keep$sqrtdetinvAri

        if (predict) {
          betahat[, i+1] <- betahat[, i]
          Yhat[, i+1] <- Yhat[, i]
          MSE[i+1] <- MSE[i]
          nselect[i+1] <- nselect[i]
        }
      }


    }
    if (display) {
      if (k %% 500 == 0) {
        print(paste("iteration is ", k))
      }
      k <- k + 1
      # print(k)
    }



  }
  if (display) {
    print("Gibbs Sampling is ending and will print the frequency of select variables.")
    #print(colSums(gamma[nburnin:(nburnin + niter), ]))
    gammaresult <- colSums(gamma[nburnin:(nburnin + niter), ])
  }
  return(list(gamma = gamma, betahat = betahat, gammaresult = gammaresult, MSE = MSE, nselect = nselect, Yhat=Yhat))
}
