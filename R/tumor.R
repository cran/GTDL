#' Tumor data
#'
#' @description Times (in days) of patients in ovarian cancer study
#' 
#' @docType data
#' @name tumor
#' @usage data(tumor)
#' @format This data frame contains the following columns:
#' \itemize{
#' \item {time:} {survival time in days}
#' \item {censured:} {censored = 0, dead = 1}
#' \item {group:} {large tumor = 0, small tumor = 1}
#' }
#' @references
#' 
#'\itemize{
#'\item Colosimo, E. A and  Giolo, S. R. Análise de Sobrevivência Aplicada. 
#'Edgard Blucher: São Paulo. 2006.
#'}
#'
#' @examples
#'
#'
#' data(tumor)
#' head(tumor)
#'
NULL