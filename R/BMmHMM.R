#' Baum-Welch Function for multivariate data
#'
#' This function is based on an EM algorithm (the Baum-Welch algorithm) and computes the parameters of a HMM for functional data, returning an object of \code{S3} class \code{mhmm}.
#' @param hmm a hmm object obtained from set.mhmm
#' @return The function returns all the parameters of a HMM, computed via Baum-Welch algorithm: the vector of initial probabilities, the transition matrix, the parameters of the distributions related to the states and the value of the log-likelihood.
#' @references
#' Martino A., Guatteri, G. and Paganoni A. M., Multivariate Hidden Markov Models for disease progression, Mox Report 59/2018, 2018
#' @import
#' graphics
#' gmfd
#' @export
#' @examples
#' data(copulahmmdata)
#' Obs <- copulahmmdata
#' n <- 20 #number of observations per statistical unit
#' n_tot <- dim(Obs)[1]
#' bt <- seq(1, n_tot, by = n)
#' distr <- c("exp", "gaussian")
#' #Initialize the HMM
#' hmm <- set_mhmm(Obs, bT = bt, nStates = 2, distr = distr)
#'
#' # Compute the parameters of the HMM with the Baum-Welch algorithm
#' bw <- fitBM_mhmm(hmm)


fitBM_mhmm <- function(hmm) {
  #retrieve initialization from default model
  tol <- 1e-2
  Obs <- hmm$Obs
  nu <- hmm$nu
  A <- hmm$A
  B <- hmm$B
  params <- hmm$params
  corr <- hmm$corr
  distr <- hmm$distr
  nStates <- length(nu)
  nObs <- dim(Obs)[1]
  bT <- hmm$bT
  eT <- c(bT[-1] - 1,nObs)
  D <- ncol(Obs)
  logLike <- 0
  logLikeOld <- 1

  while(abs(logLike-logLikeOld) > tol) {

    logLikeOld <- logLike

    #Matrix B computation
    for (i in 1:nStates) {
      B[,i] <- dmdistr(Obs, lapply(params, "[", i, TRUE), corr=corr[,,i], distr = distr)
    }
    B[which(B==0)] <- 2e-16
    B[which(is.na(B))] <- 2e-16

    #Forward-Backward application
    fb <- forwardbackward( nStates, nu, A, B, bT )
    alpha<-fb$alpha
    beta<- fb$beta
    C <- fb$c

    #Log-Likelihood computation
    logLike <- rep(0, length(bT))
    for (i in 1:length(logLike)) {
      logLike[i] <- -sum(log(C[bT[i]:eT[i]]))
    }
    logLike <- sum(logLike)

    #Xi and Gamma computation
    xi <- array(0, dim = c(nStates, nStates, nObs))
    den <- rep(0,nObs)
    for (k in 1:length(bT)) {
      for(t in bT[k]:(eT[k]-1)) {
        den[t] <- sum( (alpha[t, ] %*% A) * B[t+1, ] * beta[t+1, ]  )
      }
    }
    for (k in 1:length(bT)) {
      for(t in bT[k]:(eT[k]-1)) {
        for (i in 1:nStates) {
          for (j in 1:nStates) {
            xi[i, j, t] <- (alpha[t, i] * A[i, j] * B[t+1, j] * beta[t+1, j]) / den[t]
          }
        }
      }
    }

    gamma <- matrix( 0, nrow = nObs, ncol = nStates)
    for (k in 1:length(bT)) {
      for(t in bT[k]:(eT[k]-1)) {
        for (i in 1:nStates) {
          gamma[t, i] <- sum( xi[i, ,t] )
        }
      }
    }

    ##Parameters update

    #Initial probabilities
    nu <- colMeans(gamma[bT, ])

    #Transition matrix A
    for (i in 1:nStates) {
      for (j in 1:nStates) {
        A[i, j] <- sum(xi[i, j, ])/ sum(gamma[, i])
      }
    }

    #State-dependent distribution parameters
    temp <- array(rep(0,D*D*nObs), dim= c(D,D,nObs))
    for (d in 1:D)  {
      for (j in 1:nStates) {
        if (distr[d] == "exp") {
          params[[d]][j] <- sum( gamma[, j] ) / sum( gamma[, j] %*% Obs[,d] )
        }
        if (distr[d] == "gamma") {
          params[[d]][j, 1] <- (sum( gamma[, j] %*% Obs[,d] ) /  sum( gamma[, j] ))^2 / (sum( gamma[, j] %*% Obs[,d]^2 )/ sum( gamma[, j] ) - (sum( gamma[, j] %*% Obs[,d] ) /  sum( gamma[, j] ))^2)
          params[[d]][j, 2] <- (sum( gamma[, j] %*% Obs[,d] ) /  sum( gamma[, j] )) / (sum( gamma[, j] %*% Obs[,d]^2 )/ sum( gamma[, j] ) - (sum( gamma[, j] %*% Obs[,d] ) /  sum( gamma[, j] ))^2)
        }
        if (distr[d] == "gaussian") {
          params[[d]][j, 1] <- sum( gamma[, j] %*% Obs[,d] )/ sum( gamma[, j] )
          params[[d]][j, 2] <- sqrt( sum( gamma[, j] %*% Obs[,d]^2 )/ sum( gamma[, j] ) - (sum( gamma[, j] %*% Obs[,d] ) /  sum( gamma[, j] ))^2 )
        }
      }
    }

    #Correlation matrices
    for (j in 1:nStates) {	### Poi allo stesso modo, con le formule del Rubiner, aggiorno le stime dei parametri
      for (t in 1:nObs) {
        tempmean <- rep(0, D)
        for (d in 1:D) {
          if (distr[d] == "exp") {
            tempmean[d] <- 1/params[[d]][j]
          }
          else if (distr[d] == "gamma") {
            tempmean[d] <- params[[d]][j, 2]/params[[d]][j, 1]
          }
          else if (distr[d] == "gaussian") {
            tempmean[d] <- params[[d]][j, 1]
          }
        }
        temp[,,t] <- (Obs[t,] - tempmean)%*%t(Obs[t,] - tempmean)
      }
      corr[,,j] <- cov2cor(( rowSums((rep(gamma[,j], each = D*D) * temp ), dims = 2) / sum( gamma[, j] ) ))
    }
  }

  ## Return and class definition
  bw <- list(nu = nu, A = A, params = params, corr = corr, logLike = logLike, distr = distr, B = B, bT = bT)
  class(bw) <- "mhmm"
  return(bw)
}

