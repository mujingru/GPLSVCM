#' Fit the model with K fold cross validation and compute the CV residuals
#'
#' This is an internal function of package \code{GPLSVCM}.
#'
#' @importFrom BPST basis
#'
CV_fit=function(Y_train,X_l_train,X_nl_train,ind_l,ind_nl,U_train,X_l_pred,X_nl_pred,U_pred,B_pred,V,Tr,d,r,lambda, family, off, r_theta, eps, nfold){
  linkinv = family$linkinv
  n = length(Y_train)
  n_pred=dim(U_pred)[1]
  sfold = round(n/nfold)
  shuffle = sample(1:n)
  cv_record=rep(0,n)
  ypred=matrix(,nrow=n,ncol=n_pred)

  for(ii in 1:nfold){
    if(ii < nfold){
      vts = sort(shuffle[((ii - 1) * sfold + 1):(ii * sfold)])
    }
    if(ii == nfold){
      vts = sort(shuffle[((ii - 1) * sfold + 1):n])
    }
    vtr = setdiff(1:n, vts)

    Y_ts = as.matrix(Y_train[vts])
    Y_tr =as.matrix(Y_train[vtr])
    X_l_ts = as.matrix(X_l_train[vts, ])
    X_l_tr = as.matrix(X_l_train[vtr, ])
    X_nl_ts = as.matrix(X_nl_train[vts, ])
    X_nl_tr = as.matrix(X_nl_train[vtr, ])
    U_ts = as.matrix(U_train[vts, ])
    U_tr= as.matrix(U_train[vtr, ])
    n_ts=length(vts)

    mfit_ii = gplsvcm_est(Y_tr, X_l_tr, X_nl_tr, U_tr, V, Tr, d, r, lambda,family,off, r_theta,eps)

    if (length(ind_nl)==0){
      beta_ts=as.matrix(mfit_ii$coefficients)
      alpha_ts=as.matrix(rep(0,n_ts))
    }else{
      B_ts = suppressWarnings(basis(V,Tr,d,r,U_ts))
      alpha_ts=B_ts$B%*%mfit_ii$Qtheta
    }
    if (length(ind_l)==0){
      beta_ts=c(0)
    }
    if (length(ind_nl)!=0 & length(ind_l)!=0){
      beta_ts=mfit_ii$beta_hat
    }

    ypred_ts = linkinv(apply((X_nl_ts*alpha_ts),1,sum)+X_l_ts%*%beta_ts)
    cv_record[vts]=abs(ypred_ts-Y_ts)   #the vector for saving the K fold CV residuals for each training data point

    ####evaluate the fitted model at the test data points
    if (length(ind_nl)==0){
      beta_pred=as.matrix(mfit_ii$coefficients)
      alpha_pred=as.matrix(rep(0,dim(X_l_pred)[1]))
    }else{
      alpha_pred=B_pred%*%mfit_ii$Qtheta
    }
    if (length(ind_l)==0){
      beta_pred=c(0)
    }
    if (length(ind_nl)!=0 & length(ind_l)!=0){
      beta_pred=mfit_ii$beta_hat
    }

    y_pred= linkinv(matrix(apply((X_nl_pred*alpha_pred),1,sum)+X_l_pred%*%beta_pred,nrow=1,byrow=TRUE))  #mu_hat_-K(i)(X_pred)
    mat=y_pred[rep(1,n_ts), ]
    ypred[vts,]=mat
  }
  list(residu=cv_record,ypred=ypred)
}
