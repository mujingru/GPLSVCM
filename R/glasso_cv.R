#' @importFrom MGLM kr
#' @importFrom grpreg cv.grpreg
glasso_cv=function(Y,X,Bc,family,lambda1){
  Y=as.vector(Y)
  n=length(Y)
  J = dim(Bc)[2]
  n_p=ncol(X)
  X_m=as.matrix(cbind(rep(1,n),X))
  # The design matrix
  Z = as.matrix(cbind(X,kr(X_m, Bc, byrow=TRUE)))
  group=c(1:n_p,rep((n_p+1):(2*n_p+1), each = J))
  cvfit= cv.grpreg(Z,Y, group, penalty="grLasso",family=family$family,lambda=lambda1) #10-fold CV to select lambda
  fit=cvfit$fit   #fit on the whole data
  lambda_sel=cvfit$lambda.min    #The value of lambda with the minimum cross-validation error
  beta=fit$beta[,cvfit$min]   #the estimated coefficients corresponding to the selected lambda
  names(beta)=c("intercept",group)
  coef_l=beta[2:(n_p+1)]
  ind_l0=(1:n_p)[which(coef_l==0)]
  coef_nl=beta[-(1:(n_p+1))]
  ind_nl0=unique(names(coef_nl)[which(coef_nl==0)])
  ind_nl=as.integer(setdiff(unique(names(coef_nl)),ind_nl0))-(n_p+1)
  ind_l=setdiff(0:n_p,c(ind_nl,ind_l0))
  list(ind_nl=ind_nl,ind_l=ind_l,lambda_sel1=lambda_sel)
}
