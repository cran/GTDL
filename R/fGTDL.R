#'@name fGTDL
#'@aliases fGTDL dGTDL hGTDL sGTDL rGTDL fires
#'
#'@title The GTDL distribution
#'
#'@description Density function, survival function, failure function and random 
#'generation for the GTDL distribution.
#'@param t vector of integer positive quantile. 
#'@param param parameters (alpha and gamma are scalars, lambda non-negative).
#'@param n number of observations. 
#'@param log logical; if TRUE, probabilities p are given as log(p).
#'
#'@references
#'
#'\itemize{
#' \item Mackenzie, G. (1996). Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family. Journal of the Royal Statistical Society. 
#' Series D (The Statistician). 45. 21-34. 
#'}
#'
#'@return \code{dGTDL} gives the density function, \code{hGTDL} gives the failure function, \code{sGTDL} gives the survival function and \code{rGTDL} generates random samples.
#'@return Invalid arguments will return an error message.
#'
#'@source [d-p-q-r]GTDL are calculated directly from the definitions. 
#'
#'@details
#'\itemize{
#'\item Density function
#'\deqn{f(t\mid \boldsymbol{\theta})=\lambda\left(\frac{\exp\{\alpha{t}+\boldsymbol{X}^{\top}\boldsymbol{\beta}\}}{1+\exp\{\alpha{t}+\boldsymbol{X}^{\top}\boldsymbol{\beta}\}}\right)\times\left(\frac{1+\exp\{\alpha{t}+\boldsymbol{X}^{\top}\boldsymbol{\beta}\}}{1+\exp\{\boldsymbol{X}^{\top}\boldsymbol{\beta}\}}\right)^{-\lambda/\alpha}}
#'
#'\item Survival function
#'\deqn{S(t \mid \boldsymbol{\theta})=\left(\frac{1+\exp\{\alpha{t}+\boldsymbol{X}^{\top}\boldsymbol{\beta}\}}{1+\exp\{\boldsymbol{X}^{\top}\boldsymbol{\beta}\}}\right)^{-\lambda/\alpha}}
#' 
#'\item Failure function
#'\deqn{h(t\mid\boldsymbol{\theta})=\lambda\left(\frac{\exp\{\alpha{t}+\boldsymbol{X}^{\top}\boldsymbol{\beta}\}}{1+\exp\{\alpha{t}+\boldsymbol{X}^{\top}\boldsymbol{\beta}\}}\right)}
#' 
#' }
#' 
#'@examples
#' 
#' library(GTDL)
#' t <- seq(0,20,by = 0.1)
#' lambda <- 1.00
#' alpha <- -0.05
#' gamma <- -1.00
#' param <- c(lambda,alpha,gamma)
#' y1 <- hGTDL(t,param)
#' y2 <- sGTDL(t,param)
#' y3 <- dGTDL(t,param,log = FALSE)
#' tt <- as.matrix(cbind(t,t,t))
#' yy <- as.matrix(cbind(y1,y2,y3))
#' matplot(tt,yy,type="l",xlab="time",ylab="",lty = 1:3,col=1:3,lwd=2)
#' 
#' 
#' y1 <- hGTDL(t,c(1,0.5,-1.0))
#' y2 <- hGTDL(t,c(1,0.25,-1.0))
#' y3 <- hGTDL(t,c(1,-0.25,1.0))
#' y4 <- hGTDL(t,c(1,-0.50,1.0))
#' y5 <- hGTDL(t,c(1,-0.06,-1.6))
#' tt <- as.matrix(cbind(t,t,t,t,t))
#' yy <- as.matrix(cbind(y1,y2,y3,y4,y5))
#' matplot(tt,yy,type="l",xlab="time",ylab="Hazard function",lty = 1:3,col=1:3,lwd=2)
#' 
#' 
#' 

#'@rdname fGTDL 
#'@export

dGTDL<-function(t,param,log = FALSE){
  lambda <- param[1]
  alpha <- param[2]
  gamma <- param [3]
  d1 <- ((lambda*exp(t*alpha+gamma))/(1+exp(t*alpha+gamma)))
  d2 <- ((1+exp(t*alpha+gamma))/(1+exp(gamma)))^(-lambda/alpha)
  
  if(log == FALSE){
    return(d1*d2)
  }
  else {
    return(log(d1*d2))
    }
  }

#'@rdname fGTDL 
#'@export

hGTDL <- function(t,param){
  lambda <- param[1]
  alpha <- param[2]
  gamma <- param[3]
  h1 <- (lambda*exp(t*alpha+gamma))
  h2 <- (1+exp(t*alpha+gamma))
  return(h1/h2)
}

#'@rdname fGTDL 
#'@export

sGTDL <- function(t,param){
  lambda <- param[1]
  alpha <- param[2]
  gamma <- param[3]
  return(dGTDL(t,param)/hGTDL(t,param))
}

#'@rdname fGTDL
#'@export

rGTDL <- function(n,param){
  lambda <- param[1]
  alpha <- param[2]
  gamma <- param[3]
  u <- runif(n)
  t <- (1/alpha)*(log((1+exp(gamma))*(1-u)^(-alpha/lambda)-1)-gamma)
  return(t)
}
