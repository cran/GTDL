like1 <- function(param,t){ 
  f1 <- sum(dGTDL(param = param,t = t,log = TRUE))
  return(-f1)
}

#'@title Maximum likelihood estimation
#'
#'@description Estimate of the parameters.
#'
#'@param start Initial values for the parameters to be optimized over.
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param method The method to be used.
#'
#'@references
#'
#'\itemize{
#'\item Aarset, M. V. (1987). How to Identify a Bathtub Hazard Rate. IEEE Transactions on Reliability, 36, 106–108.
#'\item Mackenzie, G. (1996) Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family. Journal of the Royal Statistical Society. 
#' Series D (The Statistician). 45. 21-34.
#'}
#'
#'@seealso 
#'\code{\link{optim}}
#'
#'@return Returns a list of summary statistics
#' of the fitted GTDL distribution.
#' 
#'@examples
#'
#' # times data (from Aarset, 1987))
#'data(artset1987)
#'mod <- mle1.GTDL(c(1,-0.05,-1),t = artset1987)
#'
#'@export
#'@import stats
mle1.GTDL <- function(start,t,method = 'BFGS'){
  op <- suppressWarnings(optim(par = start,fn = like1,method = method,t = t,hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  rownames(mTab$Coefficients) <- c("lambda","alpha","beta")
  mTab$Coefficients <- round(mTab$Coefficients,4)
  mTab$AIC <- -2*mTab$Lik+2*(2+length(mTab$Coefficients[,1]))
  mTab$BIC <- -2*mTab$Lik+(2+length(mTab$Coefficients[,1]))*log(length(t))
  return(mTab)
}

