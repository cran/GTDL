reability.GTDL <- function(t,param){
  
  lambda <- param[1]
  alpha <- param[2]
  gamma <- param[3]
  
  func1 <- 1+exp(alpha*t+gamma)
  func2 <- 1+exp(gamma)
  pot <- -(lambda/alpha)
  
  R <- (func1/func2)^pot
  
  return(R)
}

envelope2.GTDL <- function(x){
  U	         <- x
  n	         <- length(x)
  d2s 	     <- sort(U)
  xq2 	     <- qnorm(ppoints(n))
  Xsim 	     <- matrix(0, 100, n)
  for(i in 1:100){
    u2       <- rnorm(n)
    Xsim[i,] <- u2
  }
  Xsim2      <- apply(Xsim, 1, sort)
  d21        <- matrix(0, n, 1)
  d22        <- matrix(0, n, 1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,], 0.025)
    d22[i]  <- quantile(Xsim2[i,], 0.975)
  }
  d2med      <- apply(Xsim2, 1, mean)
  fy         <- range(d2s, d21, d22)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  graphics::plot(xq2, d2s, xlab = quote("Theoretical quantiles"),
       ylab = quote("NRSP residuals"), 
       pch = 20, ylim = fy,cex.axis = 1.2,
       cex.lab = 1.2, cex = 0.6, bg = 5)
  graphics::par(new = T)
  graphics::plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
  graphics::par(new = T)
  graphics::plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "", lwd=2)
  graphics::par(new = T)
  graphics::plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
}

#'@title Normally-transformed randomized survival 
#'probability residuals for the GTDL model 
#'
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param pHat Estimate of the parameters from the GTDL model.
#'@param censur Censoring status 0=censored, a=fail.
#'@param formula The structure matrix of covariates of dimension n x p.
#'
#'@references
#'
#'\itemize{
#'\item Li, L., Wu, T., e Cindy, F. (2021). Model diagnostics for censored regression via randomized
#'survival probabilities. Statistics in Medicine, 40, 1482–1497.
#'\item de Oliveira, L. E. F., dos Santos L. S., da Silva, P. H. F., Fabio, L. C.,
#' Carrasco, J. M. F.(2022).  Análise de resíduos para o modelo logístico 
#' generalizado dependente do tempo (GTDL). Submitted. 
#'}
#'
#'@return Normally-transformed randomized survival 
#'probability residuals 
#'
#'@examples
#'
#'### Example 1
#'
#'require(survival)
#'data(lung)
#'lung <- lung[-14,]
#'lung$sex <- ifelse(lung$sex==2, 1, 0)
#'lung$ph.ecog[lung$ph.ecog==3]<-2
#'t1 <- lung$time
#'formula1 <- ~lung$sex+factor(lung$ph.ecog)+lung$age
#'censur1 <- ifelse(lung$status==1,0,1)
#'start1 <- c(0.03,0.05,-1,0.7,2,-0.1)
#'fit.model1 <- mle2.GTDL(t = t1,start = start1,
#'            formula = formula1,
#'            censur = censur1)
#'r1 <- nrsp.GTDL(t = t1,formula = formula1 ,pHat = fit.model1$Coefficients[,1],
#'              censur = censur1)
#'r1
#'
#'### Example 2
#'
#'data(tumor)
#'t2 <- tumor$time
#'formula2 <- ~tumor$group
#'censur2 <- tumor$censured
#'start2 <- c(1,-0.05,1.7)
#'fit.model2 <- mle2.GTDL(t = t2,start = start2,
#'                        formula = formula2,
#'                        censur = censur2)
#'r2 <- nrsp.GTDL(t = t2,formula = formula2, pHat = fit.model2$Coefficients[,1],
#'             censur = censur2)
#'r2

#'@export
#'@import stats
#'@import graphics
#'@import survival

nrsp.GTDL <- function(t,formula,pHat,censur){
  x.aux <- model.matrix(formula)
  x <- matrix(x.aux[,-1],ncol = (ncol(x.aux)-1))
  p <- ncol(data.matrix(x))
  conf <- NULL
  for(i in 1:length(t)){
    gamma.aux <-  x[i,]%*%matrix(pHat[3:(p+2)])
    param.aux <- c(pHat[1:2],gamma.aux)
    conf[i] <- reability.GTDL(t=t[i],param=param.aux)  
  }
  
  nrsp.r <- qnorm(censur*(conf) + (1-censur)*runif(length(t))*conf)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")
  graphics::plot(nrsp.r, xlab = "Index", ylab = "NRSP residuals",
       pch = 15, main = "", cex.axis = 1.2,
       cex.lab = 1.2, cex = 0.6, bg = 5)
  envelope2.GTDL(nrsp.r)
  return(nrsp.r)
}
