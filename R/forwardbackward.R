#' Forward-Backward Function
#'
#' This function allows you to perform the forward-backward algorithm for a Hidden Markov Model
#' @param nStates an integer representing the number of state of the Hidden Markov Model.
#' @param nu a vector of the initial probabilities.
#' @param A the transition matrix.
#' @param B a matrix containing the values of the emission probabilities
#' @param bT the vector of the beginning times for the statistical units
#' @return The function returns a numeric value indicating the distance between the two curves.
#' @keywords distance
#' @references
#' Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings of the IEEE, 77(2), 257-286.
#' @import
#' graphics
#' gmfd

forwardbackward <- function(nStates, nu, A, B, bT) {

  nObs <- dim(B)[1]
  eT <- c(bT[-1] - 1,nObs)

  #forward step
  alpha <- matrix( 0, nrow = nObs, ncol = nStates )
  c <- rep(NA, nObs)
  for (i in 1:length(bT)) {
    alpha[bT[i], ] <- nu * B[bT[i], ]
    c[bT[i]] <- 1 / sum(alpha[bT[i],])
    alpha[bT[i], ] <- alpha[bT[i], ] * c[bT[i]]
    for (t in (bT[i] + 1):eT[i]) {
      for (j in 1:nStates) {
        alpha[t, j] <- ( alpha[(t - 1), ] %*% A[, j] ) * B[t, j]
      }
      c[t] <- 1 / sum(alpha[t, ])
      alpha[t, ] <- alpha[t, ] * c[t]
    }
  }

  #backward step
  beta <- matrix( 0, nrow = nObs, ncol = nStates)

  beta[eT, ] <- rep(1,nStates)
  beta[eT, ] <- beta[eT, ] * c[eT]
  for (j in length(eT):1) {
    for (t in (eT[j] - 1):bT[j]) {
      for (i in 1:nStates) {
        beta[t, i] <- sum( beta[ t + 1, ] * A[i, ] * B[ t + 1, ] )
      }
      beta[t, ] <- beta[t, ] * c[t]
    }
  }

  #return
  return(list(alpha=alpha, beta=beta, c = c))
}
