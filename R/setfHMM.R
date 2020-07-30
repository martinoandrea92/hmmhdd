#' \code{S3} Class for Hidden Markov Models with functional response.
#'
#'  This function creates an object of \code{S3} class \code{fhmm},
#'  which contains the elements of a Hidden Markov Model with
#'  multivariate response whose components can be correlated through a
#'  gaussian copula.
#'
#' @param Obs your data, it can be a matrix (multivariate data) or a \code{funData} object (functional data)
#' @param bT the vector of the beginning times for the statistical units
#' @param nStates an integer representing the number of states of the Hidden Markov Model.
#' @param nu a vector of the initial probabilities of the Hidden Markov Model.
#' @param A the transition matrix of the Hidden Markov Model.
#' @param centroids the initialization of the centers.
#' @return The function returns an object of \code{S3} class \code{fhmm} containing the initialization and the data of your HMM
#' @references
#' Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings of the IEEE, 77(2), 257-286.
#' Martino A., Guatteri, G. and Paganoni A. M., Hidden Markov Models for multivariate functional data, MOX Report 21/2019, 2019
#' @seealso \code{\link{fitBM_fhmm}}
#' @import
#' graphics
#' gmfd
#' roahd
#' @export
#' @examples
#' \donttest{
#' data(simulatedFD)
#' FD <- simulatedFD
#' n <- 20
#' n_tot <- dim(FD$data[[1]])[1]
#' bt <- seq(1, n_tot, by = n)
#' hmm <- set_fhmm(FD, nStates = 3, bT = bt)}

set_fhmm <- function(Obs, bT = NULL, nStates = NULL, nu = 1/rep(nStates, nStates),
                     A = matrix(1/nStates, nStates, nStates), centroids = NULL) {
  nObs <- dim(Obs$data[[1]])[1]
  B <- matrix(0, nObs, nStates)
  if (is.null(centroids)) {
    centroids = funData( Obs$grid, gmfd_kmeans( Obs, n.cl = nStates, metric="L2" )$centers)
  }
  for (i in 1:nObs) {
    for (j in 1:nStates) {
      B[i, j] <- 1/funDist( funData( Obs$grid, Obs$data[[1]][i, ] ),
                           funData(centroids$grid, centroids$data[[1]][j,]), metric = "L2" )^2
    }
  }
  hmm <- list(nu = nu, A = A, B = B, centroids = centroids, bT = bT, Obs = Obs)
  class(hmm) <- "fhmm"
  return(hmm)
}
