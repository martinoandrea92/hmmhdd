#' \code{S3} Class for Hidden Markov Models with multivariate response.
#'
#'  This function creates an object of \code{S3} class \code{mhmm},
#'  which contains the elements of a Hidden Markov Model with
#'  multivariate response whose components can be correlated through a
#'  gaussian copula.
#'
#' @param Obs a n x p matrix containing the data.
#' @param bT the vector of the beginning times for the statistical units.
#' @param nStates an integer representing the number of state of the Hidden Markov Model.
#' @param nu a vector of the initial probabilities of the Hidden Markov Model.
#' @param A the transition matrix of the Hidden Markov Model.
#' @param params a list containing the state-dependent parameters of the Hidden Markov Model. Each element of the list corresponds to a component of the response vector. Each element of the list is a n x d matrix, where n is the number of states of the HMM while d is the number of parameters of the specified distribution.
#' @param corr an array containing a number of matrix equal to the number of states for the gaussian copula distribution
#' @param distr a vector containing the name of the distribution for each component, "gaussian" for the normal distribution, "gamma" for the gamma distribution, "exp" for the exponential distribution
#' @return The function returns an object of \code{S3} class \code{mhmm}, containing the initialization and the data of your HMM
#' @references
#' Martino A., Guatteri, G. and Paganoni A. M., Multivariate Hidden Markov Models for disease progression, Mox Report 59/2018, 2018
#' @seealso \code{\link{fitBM_mhmm}}
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
#'
#' #Initialize the HMM
#' hmm <- set_mhmm(Obs, bT = bt, nStates = 2, distr = distr)

set_mhmm <- function( Obs, bT = NULL, nStates = NULL, nu = rep(1/nStates, nStates),
                    A = matrix(1/nStates, nStates, nStates), corr = NULL,
                    params = NULL, distr = NULL) {
  nObs <- dim(Obs)[1]
  B <- matrix( 0, nObs, nStates)
  km <- kmeans(Obs, nStates)
  if (is.null(params)) {
    for (i in 1:dim(Obs)[2]) {
      if (distr[i] == "exp") {
        params[[i]] <- as.matrix(1/km$center[,i])
      }
      if (distr[i] == "gaussian") {
        params[[i]] <- as.matrix(cbind(km$center[,i], tapply(Obs[,i], km$cluster, sd)))
      }
      if (distr[i] == "gamma") {
        params[[i]] <- as.matrix(cbind( km$center[,i]^2/tapply(Obs[,i], km$cluster, var), km$center[,i]/tapply(Obs[,i], km$cluster, var) ) )
      }
    }
  }
  if (is.null(corr)) {
    corr <- array(diag(rep(1, dim(Obs)[2])), c(dim(Obs)[2], dim(Obs)[2], nStates))
    for (k in 1:nStates) {
      for (i in 1:dim(Obs)[2]){
        for (j in 1:dim(Obs)[2]) {
          if (i != j) {
            corr[i,j,k] <- cor( Obs[ km$cluster == k, i] , Obs[ km$cluster == k, j ] )
          }
        }
      }
    }
  }
  for (i in 1:nStates) {
    B[,i] <- dmdistr(Obs, lapply(params, "[", i, TRUE), corr = corr[,,i], distr = distr)
  }
  hmm <- list(nu = nu, A = A, B = B, params = params, corr = corr, bT = bT, distr = distr, Obs = Obs)
  class(hmm) <- "mhmm"
  return(hmm)
}
