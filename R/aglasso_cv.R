#' @importFrom MGLM kr
#' @importFrom grpreg cv.grpreg
aglasso_cv=function(Y,X,Bc,family,lambda1,lambda2){
  Y=as.vector(Y)
  n=length(Y)
  J = dim(Bc)[2]
  n_p=ncol(X)
  X_m=as.matrix(cbind(rep(1,n),X))
  # The design matrix
  Z = as.matrix(cbind(X,kr(X_m, Bc, byrow=TRUE)))
  group=c(1:n_p,rep((n_p+1):(2*n_p+1), each = J))
  cvfit0= cv.grpreg(Z,Y, group, penalty="grLasso",family=family$family,lambda=lambda1) #10-fold CV to select lambda
  fit0=cvfit0$fit   #fit on the whole data
  lambda_sel1=cvfit0$lambda.min    #The value of lambda with the minimum cross-validation error for group lasso
  beta0=fit0$beta[,cvfit0$min]   #the initial estimator from group lasso
  names(beta0)=c("intercept",group)
  coef_l=beta0[2:(n_p+1)]
  ind_l0=(1:n_p)[which(coef_l==0)]
  coef_nl=beta0[-(1:(n_p+1))]
  ind_nl0=unique(names(coef_nl)[which(coef_nl==0)])
  ind_nl=as.integer(setdiff(unique(names(coef_nl)),ind_nl0))-(n_p+1)
  ind_l=setdiff(0:n_p,c(ind_nl,ind_l0))
  ind_nlg=ind_nl
  ind_lg=ind_l
  # compute weights
  g_norm=unlist(lapply(split(beta0[-1],group),function(x) sqrt(sum(x^2))))
  weight=ifelse(g_norm==0,10^10^2.2,1/g_norm*table(group))
  # perform adaptive group lasso
  cvfit=cv.grpreg(Z,Y,group, penalty="grLasso",family=family$family,group.multiplier=weight,lambda=lambda2)
  fit=cvfit$fit   #fit on the whole data for adaptive group lasso
  lambda_sel2=cvfit$lambda.min    #The value of lambda with the minimum cross-validation error for adaptive group lasso
  beta=fit$beta[,cvfit$min]   # coefficient estimator from adptive group lasso
  names(beta)=c("intercept",group)
  coef_l=beta[2:(n_p+1)]
  ind_l0=(1:n_p)[which(coef_l==0)]
  coef_nl=beta[-(1:(n_p+1))]
  ind_nl0=unique(names(coef_nl)[which(coef_nl==0)])
  ind_nl=as.integer(setdiff(unique(names(coef_nl)),ind_nl0))-(n_p+1)
  ind_l=setdiff(0:n_p,c(ind_nl,ind_l0))
  list(ind_nlg=ind_nlg,ind_lg=ind_lg,ind_nl=ind_nl,ind_l=ind_l,lambda_sel1=lambda_sel1,lambda_sel2=lambda_sel2)
}
