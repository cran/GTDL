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

envelope.GTDL <- function(x){
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
       ylab = quote("Radomized quantile residuals"), 
       pch = 20, ylim = fy,cex.axis = 1.2,
       cex.lab = 1.2, cex = 0.6, bg = 5)
  graphics::par(new = T)
  graphics::plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
  graphics::par(new = T)
  graphics::plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "", lwd=2)
  graphics::par(new = T)
  graphics::plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
}

#'@title Randomized quantile residuals for the GTDL model
#'
#'
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param pHat Estimate of the parameters from the GTDL model.
#'@param censur censoring status 0=censored, a=fail.
#'@param formula The structure matrix of covariates of dimension n x p.
#'
#'@references
#'
#'\itemize{
#'\item Dunn, P. K. e Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational
#' and Graphical Statistics, 5, 236–244.
#'\item Louzada, F., Cuminato, J. A., Rodriguez, O. M. H., Tomazella, V. L. D., Milani, E. A., 
#'Ferreira, P. H., Ramos, P. L., Bochio, G., Perissini, I. C., Junior, O. A. G., Mota, A. L., Alegr´ıa, 
#'L. F. A., Colombo, D., Oliveira, P. G. O., Santos, H. F. L., e Magalh˜aes, M. V. C. (2020). 
#'Incorporation of frailties into a non-proportional hazard regression model and its diagnostics 
#'for reliability modeling of downhole safety valves. IEEE Access, 8, 219757 – 219774.
#'\item de Oliveira, L. E. F., dos Santos L. S., da Silva, P. H. F., Fabio, L. C.,
#' Carrasco, J. M. F.(2022).  Análise de resíduos para o modelo logístico 
#' generalizado dependente do tempo (GTDL). Submitted. 
#'}
#'
#'@details The randomized quantile residual (Dunn and Smyth, 1996), 
#'which follow a standard normal distribution is used to assess 
#'departures from the GTDL model.
#'
#'@return Randomized quantile residuals
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
#'r1 <- random.quantile.GTDL(t = t1,formula = formula1 ,pHat = fit.model1$Coefficients[,1],
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
#'r2 <- random.quantile.GTDL(t = t2,formula = formula2, pHat = fit.model2$Coefficients[,1],
#'             censur = censur2)
#'r2

#'@export
#'@import stats
#'@import graphics
#'@import survival

random.quantile.GTDL <- function(t,formula,pHat,censur){
  x.aux <- model.matrix(formula)
  x <- matrix(x.aux[,-1],ncol = (ncol(x.aux)-1))
  p <- ncol(data.matrix(x))
  conf <- NULL
  for(i in 1:length(t)){
    gamma.aux <-  x[i,]%*%matrix(pHat[3:(p+2)])
    param.aux <- c(pHat[1:2],gamma.aux)
    conf[i] <- reability.GTDL(t=t[i],param=param.aux)  
  }
  
  qr <- qnorm(censur* (1 - conf) + (1-censur)*runif(length(t),1-conf))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")
  graphics::plot(qr, xlab = "Index", ylab = "Radomized quantile residuals",
       pch = 15, main = "", cex.axis = 1.2,
       cex.lab = 1.2, cex = 0.6, bg = 5)
  envelope.GTDL(qr)
  return(qr)
}
