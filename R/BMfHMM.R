#' Baum-Welch Function for functional data
#'
#' This function is based on an EM algorithm (the Baum-Welch algorithm) and computes the parameters of a HMM for functional data, returning an object of \code{S3} class \code{fhmm}.
#' @param hmm a hmm object obtained from setfHMM
#' @return The function returns all the parameters of a HMM, computed via Baum-Welch algorithm:
#' the vector of initial probabilities, the transition matrix, the parameters of the distributions
#' related to the states and the value of the log-likelihood.
#' @references
#' Martino A., Guatteri, G. and Paganoni A. M., Hidden Markov Models for multivariate functional data, MOX Report 21/2019, 2019
#' @import
#' graphics
#' gmfd
#' roahd
#' @export
#' @examples
#' \donttest{
#' data( simulatedFD )
#' FD <- simulatedFD
#' n <- 20
#' n_tot <- dim( FD$data[[1]] )[1]
#' bt <- seq( 1, n_tot, by = n )
#' hmm <- set_fhmm( FD, nStates = 4, bT = bt )
#'
#' bw <- fitBM_fhmm(hmm)}


fitBM_fhmm <- function(hmm) {
  #retrieve initialization from default model
  tol <- 1e-2
  nu <- hmm$nu
  Obs <- hmm$Obs
  A <- hmm$A
  B <- hmm$B
  centroids <- hmm$centroids
  nStates <- length(nu)
  nObs <- dim(Obs$data[[1]])[1]
  bT <- hmm$bT
  eT <- c(bT[-1] - 1,nObs)

  logLike <- 0
  logLikeOld <- 1

  while(abs(logLike-logLikeOld) > tol) {

    logLikeOld <- logLike

    #Matrix B computation
    for (i in 1:nObs) {
      for (j in 1:nStates) {
        B[i,j] <- 1/funDist( funData(Obs$grid, lapply( Obs$data, '[', i, TRUE )), funData(Obs$grid, lapply( centroids$data, '[', j, TRUE )), metric = "L2" )^2
      }
    }

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

    #Centroids computation
    for (i in 1:nStates) {
      for (j in 1:length(centroids$data)[1]) {
        centroids$data[[j]][i,] <- rowSums( rep(gamma[,i], each=dim(Obs$data[[j]])[2]) * t(Obs$data[[j]]) ) /sum( gamma[, i] )
      }
    }
  }

  ## Return and class definition
  bw <- list(nu = nu, A = A, B = B, centroids = centroids, logLike = logLike, bT = bT, Obs = Obs)
  class(bw) <- "fhmm"
  return(bw)
}
