#' Predictions for responses of new test points from a fitted generalized
#' partial linear spatially varying coefficient model.
#'
#' \code{gplsvcm_predict} is used to make predictions for the responses of
#' predicted points from a fitted \code{gplsvcm} object.
#'
#' @details The construction of the polynomial spline functions is via
#'   \code{\link[BPST]{basis}}
#'
#' @importFrom BPST basis
#'
#' @param mfit A fitted \code{gplsvcm} object returned from function
#'   \code{gplsvcm_fit} or \code{gplsvcm_fitwMIDF}.
#' @param Xpred The design matrix for prediction (including intercept if intercept exists).
#' @param Upred The matrix of cooridinates for prediction.
#'
#' @return A vector of predicted response.
#'
#' @examples
#' # See an example of gplsvcm_fitwMIDF.
#'
#' @export

gplsvcm_predict = function(mfit,Xpred,Upred){
  if(!is.matrix(Xpred)){
    warning("The explanatory variable, Xpred, is transformed to a matrix.")
    Xpred = as.matrix(Xpred)
  }
  if(!is.matrix(Upred)){
    warning("The coordinates, Upred, is transformed to a matrix.")
    Upred = as.matrix(Upred)
  }

  ind_l=mfit$ind_l; ind_nl=mfit$ind_nl;
  family = mfit$family; linkinv = family$linkinv;


  if(identical(Upred, mfit$U)){
    if(length(mfit$ind_nl)==0){
      ypred=mfit$fitted.values
    }else{
       X_nl_pred=mfit$X_nl
       X_l_pred=mfit$X_l
       beta_hat=mfit$beta_hat
       alpha_pred=mfit$alpha_hat
       ypred= linkinv(apply((X_nl_pred*alpha_pred),1,sum)+X_l_pred%*%beta_hat) }
    } else {
      if(length(mfit$ind_nl)==0){
        ypred= linkinv(Xpred%*%mfit$beta_hat)
      }else{
      V = mfit$V; Tr = mfit$Tr; d = mfit$d; r = mfit$r; Q2 = mfit$Q2;
      Basis_full = suppressWarnings(basis(V, Tr, d, r, Upred, FALSE, FALSE))
      ind_inside_pred = Basis_full$Ind.inside
      X_nl_pred=as.matrix(Xpred[ind_inside_pred,ind_nl])
      X_l_pred=as.matrix(Xpred[ind_inside_pred,ind_l])
      if(dim(X_l_pred)[2] == 0){
        n=dim(X_nl_pred)[1]
        X_l_pred=as.matrix(rep(0,n))
      }
      theta_hat=mfit$Qtheta
      beta_hat=mfit$beta_hat
      alpha_pred=Basis_full$B%*%theta_hat
      ypred= linkinv(apply((X_nl_pred*alpha_pred),1,sum)+X_l_pred%*%beta_hat)}
    }
  return(ypred)
}
