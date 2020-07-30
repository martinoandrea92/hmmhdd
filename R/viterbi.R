#' Viterbi Function
#'
#' This function computes the best sequence of states for a Hidden Markov Model
#' @param hmm a hmm object obtained from the setHMM function
#' @return The function returns a vector \code{q} of integers indicating the best state sequence for a HMM
#' @references
#' Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings of the IEEE, 77(2), 257-286.
#' @import
#' graphics
#' gmfd
#' @export
#' @examples
#' data(copulahmmdata)
#' Obs <- copulahmmdata
#' n <- 20 #number of statistical units
#' n_tot <- dim(Obs)[1]
#' bt <- seq(1, n_tot, by = n)
#' distr <- c("exp", "gaussian")
#' #Initialize the HMM
#' parameters <- list( as.matrix(c(1,0.25)), matrix(c(3, -1, 1, 1), nrow = 2))
#' corr <- array(c(1, 0.4, 0.4, 1, 1, 0.1, 0.1, 1), dim = c(2, 2, 2))
#' hmm <- set_mhmm(Obs, bT = bt, nStates = 2, params = parameters, corr = corr, distr = distr)
#' # Compute the parameters of the HMM with the Baum-Welch algorithm
#' bw <- fitBM_mhmm(hmm)
#'
#' v <- viterbi(bw)

viterbi <- function(hmm) {
  #Retrieve initialization from B-W output
  nu <- hmm$nu
  A <- hmm$A
  B <- hmm$B
  bT <- hmm$bT
  nObs <- nrow(B)
  eT <- c(bT[-1] - 1,nObs)
  nStates <- ncol(B)
  delta <- matrix(0, nrow = nObs, ncol = nStates)
  psi <- matrix(0, nrow = nObs, ncol = nStates)

  ## Initialization
  delta[1, ] <- nu * B[1, ]
  for (k in 1:length(bT)) {
    delta[bT[k], ] <- nu * B[bT[k], ]
  }

  ## Recursion step
  M <- NULL
  for (k in 1:length(bT)) {
    for(t in (bT[k]+1):(eT[k])) {
      temp <- rep(0, nStates)
      for (j in 1:nStates) {
        temp <- rep(0, nStates)
        for (i in 1:nStates) {
          temp[i] <- delta[t - 1, i] * A[i, j]
        }
        M <- max(temp)
        psi[t, j] <- which(temp == M)
        delta[t, j] <- M * B[t, j]
      }
    }
  }

  P <- rep(0, length(bT))
  q <- rep(0,nObs)
  for (k in 1:length(bT)) {
    P[k]<- max(delta[eT[k], ])
    q[eT[k]] <- which(delta[eT[k], ] == P[k])
  }

  for (k in 1:length(bT)) {
    for(t in (eT[k]-1):(bT[k])) {
      q[t] <- psi[t+1, q[t+1]]
    }
  }

  #Return
  return(q)
}
