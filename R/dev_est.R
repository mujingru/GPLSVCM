#' Calculate the generalized deviance of a fitted model
#'
#' \code{dev_est} is an internal function of \code{gplsvcm_fitwMIDF}.
#'
#' @param X_l The matrix of linear covariates for observations.
#' @param X_nl The matrix of nonlinear covariates for observations.
#' @param mfit A list of information returned by function \code{gplsvcm_est}
#' @param family The family object which specifies the distribution and link to
#'   use.
#'
dev_est=function(X_l,X_nl,mfit,family){
  if(dim(X_l)[2] == 0){
    return(dev1(mfit,family))
  }
  if(dim(X_nl)[2] == 0){
    return(dev2(X_l,mfit,family))
  }
  if(dim(X_l)[2] != 0 & dim(X_nl)[2] != 0){
    return(dev3(X_l,X_nl,mfit,family))
  }
}

