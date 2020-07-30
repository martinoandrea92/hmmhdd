

#' Copula density Function for multivariate data
#'
#' This function returns the value of a multivariate density defined on a gaussian copula
#' @param x a n x d matrix containing the multivariate data
#' @param params a list of D elements corresponding to the   matrices corresponding to the parameters:
#' @param corr a matrix containing the correlation parameters of the multivariate distribution related to the states of the HMM;
#' @param distr a vector containing the distribution for each component, "gaussian" for the normal distribution, "gamma" for the gamma distribution, "exp" for the exponential distribution
#' @return The function returns a numeric value corresponding to the value of the density of the multivariate distribution where the dependence is computed using a copula.
#' @references  Song, P. X.-K.(2000). Multivariate dispersion models generated from Gaussian copula. Scandinavian Journal of Statistics, 27(2), 305-320.
#' @import
#' stats
#' mvtnorm
#' @export
#' @examples
#' #Load the copuladata
#' data(copuladata)
#' params <- list(t(c(1.5, 1)), as.matrix(0.5))
#' corr <- cbind(c(1,0.7), c(0.7,1))
#' distr <- c("gaussian", "exp")
#'
#' density <- dmdistr(copuladata, params = params, corr = corr, distr = distr)

dmdistr <- function(x, params, corr=diag(ncol(x)), distr)
{
  # checks

  if(is.data.frame(x))
    x = as.matrix(x)
  else if (is.vector(x))
    x = matrix(x, nrow=1)

  D = ncol(x)
  cdf <- matrix(0, dim(x)[1], D)
  for (d in 1:D) {
    if (distr[d] == "exp") {
      cdf[,d] <- pexp(x[,d], rate = params[[d]])
    }
    else if (distr[d] == "gamma") {
      cdf[,d] <- pgamma(x[,d], shape = params[[d]][1], rate = params[[d]][2])
    }
    else if (distr[d] == "gaussian") {
      cdf[,d] <- pnorm(x[,d], mean = params[[d]][1], sd = params[[d]][2])
    }
  }
  normcdf = qnorm(cdf)

  ##copula contribution to log-PDF

  res.copula.add = dmvnorm(normcdf, sigma=corr, log=TRUE)
  res.copula.sub = ifelse(dim(x)[1]>1, rowSums(dnorm(normcdf, log=TRUE)), sum(dnorm(normcdf, log=TRUE)))
  res.copula = res.copula.add - res.copula.sub

  ## marginal contribution to log-PDF
  md <- matrix(0, dim(x)[1], D)
  for (d in 1:D) {
    if (distr[d] == "exp") {
      md[,d] <- dexp(x[,d], rate = params[[d]], log = T)
    }
    else if (distr[d] == "gamma") {
      md[,d] <- dgamma(x[,d], shape = params[[d]][1], rate = params[[d]][2], log = T)
    }
    else if (distr[d] == "gaussian") {
      md[,d] <- dnorm(x[,d], mean = params[[d]][1], sd = params[[d]][2], log = T)
    }
  }
  res.data = rowSums(md)

  ## return

  retn = res.copula + res.data

  return(exp(retn))
}


#' Random generator Function for multivariate data
#'
#' This function generates the values of a multivariate distribution with given marginals and correlation matrix
#' @param n the length of the sample
#' @param params a list of D elements corresponding to the   matrices corresponding to the parameters:
#' @param corr a matrix containing the correlation parameters of the multivariate distribution related to the states of the HMM;
#' @param distr a vector containing the distribution for each component, "gaussian" for the normal distribution, "gamma" for the gamma distribution, "exp" for the exponential distribution
#' @return The function returns a n x d matrix
#' @references  Song, P. X.-K.(2000). Multivariate dispersion models generated from Gaussian copula. Scandinavian Journal of Statistics, 27(2), 305-320.
#' @import
#' stats
#' mvtnorm
#' @export
#' @examples
#' #Simulate a bivariate sample with a gaussian component and an exponential one
#' n <- 100
#' params <- list(t(c(1,1)), as.matrix(1))
#' corr <- cbind(c(1,0.5), c(0.5,1))
#' distr <- c("gaussian", "exp")
#'
#' x <- rmdistr(n = n, params = params, corr = corr, distr = distr)

rmdistr <- function(n, params, corr, distr)
{
  # checks
  if((!is.matrix(corr) || !isSymmetric(corr)))
    stop("'corr' must be a symmetric matrix")
  if(!is.list(params))
    params <- list(params)
  D <- length(distr)
  # univariate case
  if(D == 1) {
    if (distr[D] == "exp") {
      return(rexp(n, rate = params[[D]]))
    }
    else if (distr[D] == "gamma") {
      return(rgamma(n, shape = params[[D]][1], rate = params[[D]][2]))
    }
    else if (distr[D] == "gaussian") {
      return(rnorm(n, mean = params[[D]][1], sd = params[[D]][2]))
    }
  }
  ## generate standard multivariate normal matrix and convert to CDF

  Z = rmvnorm(n, sigma=corr)
  cdf = pnorm(Z)

  #convert each marginal to its distribution
  md <- matrix(0, n, D)
  for (d in 1:D) {
    if (distr[d] == "exp") {
      md[,d] <- qexp(cdf[,d], rate = params[[d]])
    }
    else if (distr[d] == "gamma") {
      md[,d] <- qgamma(cdf[,d], shape = params[[d]][1], rate = params[[d]][2])
    }
    else if (distr[d] == "gaussian") {
      md[,d] <- qnorm(cdf[,d], mean = params[[d]][1], sd = params[[d]][2])
    }
  }

  #return
  return(md)
}
