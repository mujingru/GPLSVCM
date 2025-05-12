#' Fit the model with the i-th training data point removed and compute the leave
#' one out residuals
#'
#' This is an internal function of package \code{GPLSVCM}.
#'
#' @importFrom BPST basis

loo_fit=function(Y_train,X_l_train,X_nl_train,ind_l,ind_nl,U_train,X_l_pred,X_nl_pred,U_pred,B_pred,V,Tr,d,r,lambda, family, off, r_theta, eps){
  linkinv = family$linkinv
  n=length(Y_train)
  n_pred=dim(U_pred)[1]
  res=c()
  ypred=matrix(,nrow=n,ncol=n_pred)

  for (i in 1:n){
    Y_i=Y_train[-i]
    X_nli=X_nl_train[-i,]
    X_li=X_l_train[-i,]
    U_i=U_train[-i,]
    mfit_i=gplsvcm_est(Y_i,X_li,X_nli,U_i,V,Tr,d,r,lambda,family,off, r_theta,eps) #gaussian

    if (length(ind_nl)==0){
      beta_hati=as.matrix(mfit_i$coefficients)
      alpha_hati=as.matrix(c(0))
    }else{
      B_predi = suppressWarnings(basis(V,Tr,d,r,U_train[i,],FALSE,FALSE))
      alpha_hati=B_predi$B%*%mfit_i$Qtheta
    }
    if (length(ind_l)==0){
      beta_hati=c(0)
    }
    if (length(ind_nl)!=0 & length(ind_l)!=0){
      beta_hati=mfit_i$beta_hat
    }

    ypred_i= linkinv(sum(X_nl_train[i,]*alpha_hati)+X_l_train[i,]%*%beta_hati)
    residual_i=abs(Y_train[i]-ypred_i)         # the i-th leave-one-out residual
    res=c(res,residual_i)       # the vector for saving the leave-one-out residuals

    ####evaluate the fitted model at the test data points
    if (length(ind_nl)==0){
      beta_pred=as.matrix(mfit_i$coefficients)
      alpha_pred=as.matrix(rep(0,dim(X_l_pred)[1]))
    }else{
      alpha_pred=B_pred%*%mfit_i$Qtheta
    }
    if (length(ind_l)==0){
      beta_pred=c(0)
    }
    if (length(ind_nl)!=0 & length(ind_l)!=0){
      beta_pred=mfit_i$beta_hat
    }

    y_pred= linkinv(apply((X_nl_pred*alpha_pred),1,sum)+X_l_pred%*%beta_pred)  #mu_hat_-i(X_pred)
    ypred[i,]=y_pred
  }

  list(residu=res,ypred=ypred)
}
