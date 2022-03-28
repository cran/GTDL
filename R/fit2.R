like2 <- function(t,formula,censur,para){
  x_aux <- model.matrix(formula)
  x <- matrix(x_aux[,-1],ncol = (ncol(x_aux)-1))
  p <- ncol(data.matrix(x))
  ll <- NULL
  for(i in 1:dim(x)[1]){
    gama_aux <- x[i,]%*%matrix(para[3:(p+2)])
    para_aux <- c(para[1:2],gama_aux)  
    l <- (hGTDL(t = t[i],param = para_aux)^censur[i])*sGTDL(t = t[i],param = para_aux)
    ll[i] <- log(l)
  }
  return(-sum(ll))  
}

#'@title Maximum likelihood estimates of the GTDL model
#'
#'@param start Initial values for the parameters to be optimized over.
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param censur  censoring status 0=censored, a=fail.
#'@param formula The structure matrix of covariates of dimension n x p.
#'@param method The method to be used.
#'
#'@references
#'
#'\itemize{
#' \item Mackenzie, G. (1996) Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family. Journal of the Royal Statistical Society. 
#' Series D (The Statistician). (45). 21-34.
#' }
#' 
#'@seealso 
#'\code{\link{optim}}
#'
#'@return Returns a list of summary statistics
#' of the fitted GTDL model.
#'
#'@examples
#'
#'### Example 1
#'
#'require(survival)
#'data(lung)
#'
#'lung <- lung[-14,]
#'lung$sex <- ifelse(lung$sex==2, 1, 0)
#'lung$ph.ecog[lung$ph.ecog==3]<-2
#'t1 <- lung$time
#'start1 <- c(0.03,0.05,-1,0.7,2,-0.1)
#'formula1 <- ~lung$sex+factor(lung$ph.ecog)+lung$age
#'censur1 <- ifelse(lung$status==1,0,1)
#'fit.model1 <- mle2.GTDL(t = t1,start = start1,
#'                      formula = formula1,
#'                      censur = censur1)
#'fit.model1
#'
#'### Example 2
#'
#'data(tumor)
#'t2 <- tumor$time
#'start2 <- c(1,-0.05,1.7)
#'formula2 <- ~tumor$group
#'censur2 <- tumor$censured
#'fit.model2 <- mle2.GTDL(t = t2,start = start2,
#'                        formula = formula2,
#'                        censur = censur2)
#'fit.model2
#'
#'@export
#'@import stats
#'@import survival
mle2.GTDL <- function(t,start,formula,censur,method = "BFGS"){
   op <- suppressWarnings(optim(par = start,fn = like2,
                               t = t,
                               method = method,
                               formula = formula,
                               censur = censur,
                               hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  rownames(mTab$Coefficients) <- c("lambda","alpha",paste0("beta ",c(1:(dim(mTab$Coefficients)[1]-2))))
  mTab$Lik <- round(mTab$Lik,4)
  mTab$Coefficients <- round(mTab$Coefficients,4)
  mTab$AIC <- -2*mTab$Lik+2*(2+length(mTab$Coefficients[,1]))
  mTab$BIC <- -2*mTab$Lik+(2+length(mTab$Coefficients[,1]))*log(length(t))
  return(mTab)
  
}
