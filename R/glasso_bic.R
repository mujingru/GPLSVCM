#' @importFrom MGLM kr
#' @importFrom grpreg grpreg
glasso_bic=function(Y,X,Bc,family,lambda1){
  Y=as.vector(Y)
  n=length(Y)
  J = dim(Bc)[2]
  n_p=ncol(X)
  X_m=as.matrix(cbind(rep(1,n),X))
  # The design matrix
  Z = as.matrix(cbind(X,kr(X_m, Bc, byrow=TRUE)))
  group=c(1:n_p,rep((n_p+1):(2*n_p+1), each = J))
  fit0= grpreg(Z,Y,group,penalty="grLasso",family=family$family,lambda=lambda1)
  df0=colSums(as.matrix(fit0$beta[-1,])!=0)
  BIC0=log(fit0$loss)+log(n)*df0/n
  beta0=fit0$beta[,which.min(BIC0)] #the initial estimator from group lasso
  lambda_sel1=lambda1[which.min(BIC0)]
  names(beta0)=c("intercept",group)
  coef_l=beta0[2:(n_p+1)]
  ind_l0=(1:n_p)[which(coef_l==0)]
  coef_nl=beta0[-(1:(n_p+1))]
  ind_nl0=unique(names(coef_nl)[which(coef_nl==0)])
  ind_nl=as.integer(setdiff(unique(names(coef_nl)),ind_nl0))-(n_p+1)
  ind_l=setdiff(0:n_p,c(ind_nl,ind_l0))
  list(ind_nl=ind_nl,ind_l=ind_l,lambda_sel1=lambda_sel1)
}
