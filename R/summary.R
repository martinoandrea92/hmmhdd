#' Summarizing multivariate Hidden Markov Models
#'
#'
#' \code{summary} method for class "\code{mhmm}".
#' @param object an object of class "\code{mhmm}", a result of a call to \code{fitBM_mhmm}
#' @param ... additional arguments affecting the summary produced.
#' @return The function "\code{summary.mhmm}" returns a list of summary statistics for a hidden markov model obtained through the Baum-Welch algorithm.
#' @method summary mhmm
#' @references
#' Martino A., Guatteri, G. and Paganoni A. M., Multivariate Hidden Markov Models for disease progression, Mox Report 59/2018, 2018
#' @seealso \code{\link{fitBM_mhmm}}
#' @import
#' stats
#' mvtnorm
#' @export
#' @examples
#' data(copulahmmdata)
#' Obs <- copulahmmdata
#' n <- 20 #number of statistical units
#' n_tot <- 5000 #total number of observations
#' bt <- seq(1, n_tot, by = n)
#' distr <- c("exp", "gaussian")
#' #Initialize the HMM
#' parameters <- list(as.matrix(c(1,0.25)), matrix(c(3, -1, 1, 1), nrow = 2))
#' corr <- array(c(1, 0.4, 0.4, 1, 1, 0.1, 0.1, 1), dim = c(2, 2, 2))
#' hmm <- set_mhmm(Obs, bT = bt, nStates = 2, params = parameters, corr = corr, distr = distr)
#' # Compute the parameters of the HMM with the Baum-Welch algorithm
#' bw <- fitBM_mhmm(hmm)
#'
#' summary(bw)

summary.mhmm <- function( object, ... ) {
  ###Initial state probabilites###
  states <- "S_"
  output = rep( NULL, length( object$nu ) )
  for( i in 1:length( object$nu ) ) {
    Names <- sapply( states, function( s ) { paste( s, i, sep="") } )
    output <- rbind( output, Names )
  }
  vec <- as.array( object$nu )
  dimnames( vec ) <- list( output )
  cat( "Initial state probabilities:\n" )
  print( vec )
  ###Transition Matrix###

  cat( "Transition matrix:\n" )
  mat <- object$A
  rows <- "fromS_"
  cols <- "toS_"
  output1 = rep( NULL,length( object$nu ) )
  output2 = rep( NULL,length( object$nu ) )
  for( i in 1:length( object$nu ) ) {
    rowNames <- sapply( rows, function( s ) { paste( s, i, sep = "" ) } )
    colNames <- sapply( cols, function( s ) { paste( s, i, sep = "" ) } )
    output1 <- rbind( output1, rowNames )
    output2 <- rbind( output2, colNames )
  }
  dimnames( mat ) <- list( output1, output2 )
  print( mat )

  ### Response parameters ###

  for ( i in 1:length( object$distr ) ) {
    cat( "Response", i, ":" )
    if( object$distr[ i ] == "exp" )
      cat( "Exponential\n" )
    if( object$distr[ i ] == "gaussian" )
      cat( "Gaussian\n" )
    if( object$distr[i] == "gamma" )
      cat( "Gamma\n" )
  }
  cat( "Response parameters\n\n" )
  respmat <- rep( NULL, D )
  paramnames <- NULL
  gaussiannames <- c( "mean.Comp_", "sd.Comp_" )
  expnames <- c( "rate.Comp_" )
  gammanames <- c( "rate.Comp_", "shape.Comp_" )
  for ( i in 1:length( object$distr ) ) {
    if( object$distr[ i ] == "exp" ) {
      respmat <- cbind( respmat, object$params[[ i ]] )
      expNames <- sapply( expnames, function( s ) { paste( s, i, sep="" ) } )
      paramnames <- c( paramnames, expNames )
    }

    if( object$distr[ i ] == "gaussian" ) {
      respmat <- cbind( respmat, object$params[[ i ]][ ,1 ], object$params[[ i ]][, 2 ] )
      gaussianNames <- sapply( gaussiannames, function( s ) { paste( s, i, sep = "") } )
      paramnames <- c( paramnames, gaussianNames )
    }

    if( object$distr[ i ] == "gamma" ) {
      respmat <- cbind( respmat, object$params[[ i ]][ ,1 ], object$params[[ i ]][ ,2 ] )
      gammaNames <- sapply( gammanames, function( s ) { paste( s, i, sep = "") } )
      paramnames <- c( paramnames, gammaNames )
    }
  }
  colnames( respmat ) <- paramnames
  rownames( respmat ) <- output
  print( respmat )

  for ( i in 1:length( object$nu ) ) {
    cat( "Correlation matrix for State", i, "\n" )
    print( object$corr[ , ,i ] )
  }
}


#' Summarizing functional Hidden Markov Models
#'
#'
#' \code{summary} method for class "\code{fhmm}".
#' @param object an object of class "\code{fhmm}", a result of a call to \code{fitBM_fhmm}
#' @param ... additional arguments affecting the summary produced.
#' @return The function "\code{summary_fhmm}" returns a list of summary statistics for a hidden markov model with functional response obtained through the Baum-Welch algorithm.
#' @method summary fhmm
#' @references
#' Martino A., Guatteri, G. and Paganoni A. M., Hidden Markov Models for multivariate functional data, MOX Report 21/2019, 2019
#' @seealso \code{\link{fitBM_fhmm}}
#' @import
#' stats
#' @export
#' @examples
#' \donttest{
#' data(simulatedFD)
#' n <- 20
#' n_tot <- 2000
#' bt <- seq(1,n_tot,by=n)
#' FD <- simulatedFD
#' hmm <- set_fhmm(FD, nStates = 3, bT=bt)
#'
#' bw <- fitBM_fhmm(hmm)
#' summary(bw)}

summary.fhmm <- function( object, ... ) {
  ###Initial state probabilites###
  states <- "S_"
  output <- rep( NULL, length( object$nu ) )
  for(i in 1:length( object$nu ) ) {
    Names <- sapply( states, function( s ) { paste( s, i, sep = "" ) } )
    output <- rbind( output, Names )
  }
  vec <- as.array( object$nu )
  dimnames( vec ) <- list( output )
  cat( "Initial state probabilities:\n" )
  print( vec )
  ###Transition Matrix###

  cat( "Transition matrix:\n" )
  mat <- object$A
  rows <- "fromS_"
  cols <- "toS_"
  output1 <- rep( NULL, length( object$nu ) )
  output2 <- rep( NULL, length( object$nu ) )
  for(i in 1:length( object$nu ) ) {
    rowNames <- sapply( rows, function( s ) { paste( s, i, sep = "" ) } )
    colNames <- sapply( cols, function( s ) { paste( s, i, sep="" ) } )
    output1 <- rbind( output1, rowNames )
    output2 <- rbind( output2, colNames )
  }
  dimnames( mat ) <- list( output1, output2 )
  print( mat )
  cat( "\n AIC:\n" )
  print( aic( object ) )
  cat( "\n BIC:\n" )
  print( bic( object) )
}


#' AIC
#' Aikake Information Criterion for the HMM
#' @param hmm an object of class "\code{mhmm}" or "\code{fhmm}", a result of a call to \code{\link{fitBM_mhmm}} or \code{\link{fitBM_fhmm}}
#' @return The function returns the value corresponding to the AIC for a HMM.
#' @references
#' Martino A., Guatteri, G. and Paganoni A. M., Multivariate Hidden Markov Models for disease progression, Mox Report 59/2018, 2018
#' Martino A., Guatteri, G. and Paganoni A. M., Hidden Markov Models for multivariate functional data, MOX Report 21/2019, 2019
#' @seealso \code{\link{bic}}, \code{\link{fitBM_fhmm}}, \code{\link{fitBM_mhmm}}
#' @import
#' stats
#' mvtnorm
#' @export
#' @examples
#' data(copulahmmdata)
#' Obs <- copulahmmdata
#' n <- 20 #number of statistical units
#' n_tot <- 5000 #total number of observations
#' bt <- seq(1, n_tot, by = n)
#' distr <- c("exp", "gaussian")
#' #Initialize the HMM
#' parameters <- list(as.matrix(c(1,0.25)), matrix(c(3, -1, 1, 1), nrow = 2))
#' corr <- array(c(1, 0.4, 0.4, 1, 1, 0.1, 0.1, 1), dim = c(2, 2, 2))
#' hmm <- set_mhmm(Obs, bT = bt, nStates = 2, params = parameters, corr = corr, distr = distr)
#' # Compute the parameters of the HMM with the Baum-Welch algorithm
#' bw <- fitBM_mhmm(hmm)
#'
#' AIC <- aic(bw)

aic <- function(hmm) {
  logLike <- hmm$logLike
  if (class(hmm) == "mhmm") {
    cont <- 0
    for (i in 1:length(hmm$distr)) {
      cont <- cont + prod(dim(hmm$params[[i]]))
    }
    npar <- (length(hmm$nu) - 1) + dim(hmm$A)[1] * (dim(hmm$A)[1] - 1) + cont
  }
  else {
    npar <- (length(hmm$nu) - 1) + dim(hmm$A)[1] * (dim(hmm$A)[1] - 1) + length(hmm$centroids)
  }
  return(-2*logLike +  2*npar)
}


#' BIC
#' Bayesian Information Criterion for the HMM
#' @param hmm an object of class "\code{mhmm}" or "\code{fhmm}", a result of a call to \code{\link{fitBM_mhmm}} or \code{\link{fitBM_fhmm}}
#' @return The function returns the value corresponding to the BIC for a HMM.
#' @references
#' Martino A., Guatteri, G. and Paganoni A. M., Multivariate Hidden Markov Models for disease progression, Mox Report 59/2018, 2018
#' Martino A., Guatteri, G. and Paganoni A. M., Hidden Markov Models for multivariate functional data, MOX Report 21/2019, 2019
#' @seealso \code{\link{aic}}, \code{\link{fitBM_fhmm}}, \code{\link{fitBM_mhmm}}
#' @import
#' stats
#' mvtnorm
#' @export
#' @examples
#' data(copulahmmdata)
#' Obs <- copulahmmdata
#' n <- 20 #number of statistical units
#' n_tot <- 5000 #total number of observations
#' bt <- seq(1, n_tot, by = n)
#' distr <- c("exp", "gaussian")
#' #Initialize the HMM
#' parameters <- list(as.matrix(c(1,0.25)), matrix(c(3, -1, 1, 1), nrow = 2))
#' corr <- array(c(1, 0.4, 0.4, 1, 1, 0.1, 0.1, 1), dim = c(2, 2, 2))
#' hmm <- set_mhmm(Obs, bT = bt, nStates = 2, params = parameters, corr = corr, distr = distr)
#' # Compute the parameters of the HMM with the Baum-Welch algorithm
#' bw <- fitBM_mhmm(hmm)
#'
#' BIC <- bic(bw)

bic <- function(hmm) {
  logLike <- hmm$logLike
  if (class(hmm) == "mhmm") {
    cont <- 0
    for (i in 1:length(hmm$distr)) {
      cont <- cont + prod(dim(hmm$params[[i]]))
    }
    npar <- (length(hmm$nu) - 1) + dim(hmm$A)[1] * (dim(hmm$A)[1] - 1) + cont
  }
  else {
    npar <- (length(hmm$nu) - 1) + dim(hmm$A)[1] * (dim(hmm$A)[1] - 1) + length(hmm$centroids)
  }
  return(-2*logLike+log(dim(hmm$B)[1])*npar)
}



#' Plotting functional Hidden Markov Models
#'
#'
#' \code{plot} method for class "\code{fhmm}".
#' @param x an object of class "\code{fhmm}", a result of a call to \code{fitBM_fhmm}
#' @param ... additional arguments affecting the summary produced.
#' @return The function returns a plot of the centers representing the states of the Hidden Markov Model.
#' @method plot fhmm
#' @references
#' Martino A., Guatteri, G. and Paganoni A. M., Hidden Markov Models for multivariate functional data, MOX Report 21/2019, 2019
#' @import
#' stats
#' @seealso \code{\link{fitBM_fhmm}}
#' @export
#' @examples
#' \donttest{
#' data(simulatedFD)
#' n <- 20
#' n_tot <- 2000
#' bt <- seq(1,n_tot,by=n)
#' FD <- simulatedFD
#' hmm <- set_fhmm(FD, nStates = 3, bT=bt)
#' bw <- fitBM_fhmm(hmm)
#'
#' plot(bw)}

plot.fhmm <- function(x, ...) {
  C <- x$centroids$data
  R <- length(C)
  opar <- par(mfrow = c(1, R))
  on.exit(par(opar))
  grid <- x$Obs$grid
  data <- x$Obs$data
  nStates <- length(x$nu)
  for (r in 1:R) {
    matplot(grid, t(data[[r]]), type = "l", col = "grey",
            main = paste("State centers for component ", r),
            lty = 3, xlab = "t", ylab = "X(t)")
    for (k in 1:nStates) {
      lines(grid, C[[r]][k, ], col = 1 + k, lwd = 3)
    }
  }
}

