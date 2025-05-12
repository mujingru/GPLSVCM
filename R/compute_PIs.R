#' Compute the prediction intervals for responses of new test points from a
#' fitted generalized partial linear spatially varying coefficient model.
#'
#' \code{compute_PIs} compute the prediction intervals for responses of new test
#' points from a fitted \code{gplsvcm} object based on a selected prediction
#' method among Jackknife, Jackknife+ and K-fold cross validation (CV+), and
#' return the prediction intervals for responses of predicted points.
#'
#' @details The construction of the polynomial spline functions is via
#'   \code{\link[BPST]{basis}}.
#'
#' @importFrom MGLM kr
#' @importFrom BPST basis
#' @importFrom MASS glm.nb
#' @importFrom Matrix bdiag
#'
#' @param Y_train The response variable,a \code{n} by one matrix where \code{n}
#'   is the number of observations in the training data set .
#' @param X_train The design matrix of \code{n} by \code{np} where \code{np} is
#'   the number of covariates (including intercept if intercept exists). Each row is a vector of the covariates for an
#'   observation in the training data set.
#' @param ind_l The vector of the indexes which indicate the columns of linear
#'   covariates in \code{X_train}.
#' @param ind_nl The vector of the indexes which indicate the columns of
#'   nonlinear covariates in \code{X_train}.
#' @param U_train A \code{n} by two matrix where each row is the coordinates of
#'   an observation in the training data set.
#' @param X_pred The design matrix for prediction (including intercept if intercept exists).
#' @param U_pred The matrix of coordinates for prediction.
#' @param V A \code{N} by two matrix of vertices of a triangulation, where
#'   \code{N} is the number of vertices and each row is the coordinates for a
#'   vertex.
#' @param Tr A \code{n_Tr} by three triangulation matrix, where \code{n_Tr} is
#'   the number of triangles in the triangulation and each row is the indices of
#'   vertices in \code{V}.
#' @param d The degree of piecewise polynomials -- default is 2.
#' @param r The smoothness parameter and \code{r} \eqn{<} \code{d} -- default is
#'   1.
#' @param lambda The vector of the candidates of penalty parameter -- default is
#'   grid points of 10 to the power of a sequence from -6 to 6 by 0.5.
#' @param family The family object which specifies the distribution and link to
#'   use (see \code{\link{glm}} and \code{\link{family}}).
#' @param off The offset -- default is 0.
#' @param r_theta The vector of the upper and lower bound of an interval to
#'   search for an additional parameter \code{theta} for negative binomial
#'   scenario -- default is c(2,8).
#' @param eps The error tolerance for the Pearson estimate of the scale
#'   parameter, which is as close  to 1, when estimating an additional parameter
#'   \code{theta} for negative binomial scenario -- default is 0.01.
#' @param method The prediction method used in the computation, options are
#'   "CV+", "Jackknife" and "Jackknife+" -- default is "CV+".
#' @param cp The desired coverage level for the prediction intervals -- default
#'   is 0.95.
#' @param  nfold The number of folds for CV+ method -- default is 10.
#'
#' @return  A data frame of computed prediction intervals for responses of the
#'   predicted points.
#'
#' @references Barber et al.(2021) Predictive inference with the jackknife+.
#'   Ann.Statist.49(1):486-507
#'   \url{https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-1/Predictive-inference-with-the-jackknife/10.1214/20-AOS1965.full}
#'
#' @examples
#' # See an example of gplsvcm_fitwMIDF.
#'
#' @export

compute_PIs=function(Y_train, X_train,ind_l,ind_nl,U_train,X_pred,U_pred, V, Tr, d = 2, r = 1,
                     lambda = 10^seq(-6, 6, by = 0.5),family, off = 0,
                     r_theta = c(2, 8), eps= 0.01,method="CV+", cp=0.95, nfold = 10)
{
  if(!is.matrix(X_pred)){
    warning("The explanatory variable, X_pred, is transformed to a matrix.")
    X_pred = as.matrix(X_pred)
  }
  if(!is.matrix(U_pred)){
    warning("The coordinates,U_pred, is transformed to a matrix.")
    U_pred = as.matrix(U_pred)
  }
  ###preparation for fitting and evaluating the model
  Basis_full = suppressWarnings(basis(V, Tr, d, r, U_train,FALSE,FALSE))
  ind_inside = Basis_full$Ind.inside
  Y_train = as.matrix(Y_train[ind_inside])
  X_nl_train = as.matrix(X_train[ind_inside,ind_nl])
  X_l_train=as.matrix(X_train[ind_inside, ind_l])
  U_train = as.matrix(U_train[ind_inside, ])

  if (length(ind_nl)!=0){
  Ball_pred = suppressWarnings(basis(V,Tr,d,r,U_pred,FALSE,FALSE))
  ind_inside = Ball_pred$Ind
  X_nl_pred = as.matrix(X_pred[ind_inside,ind_nl])
  if (length(ind_l)==0){
  X_l_pred=as.matrix(rep(0,dim(X_nl_pred)[1]))
  }else{
  X_l_pred=as.matrix(X_pred[ind_inside, ind_l])}
  U_pred = as.matrix(U_pred[ind_inside, ])
  B_pred = as.matrix(Ball_pred$Bi)
  }else{
    X_l_pred = as.matrix(X_pred)
    X_nl_pred= as.matrix(rep(0,dim(X_l_pred)[1]))
  }
  ###compute the prediction intervals for test data points with the specified method
  if (method == "CV+"){
    cvfit=CV_fit(Y_train,X_l_train,X_nl_train,ind_l,ind_nl,U_train,X_l_pred,X_nl_pred,U_pred,B_pred,V,Tr,d,r,lambda,family,off, r_theta,eps,nfold)
    mat1=as.matrix(cvfit$residu)
    mat2=mat1[, rep(1, dim(cvfit$ypred)[2])]
    mat3=cvfit$ypred-mat2
    mat4=cvfit$ypred+mat2
    LL= apply(mat3,2,function(x) quantile(x,probs = 1-cp,type=3))
    UU= apply(mat4,2,function(x) quantile(x,probs = cp,type=3))
  }else{

    lfit=loo_fit(Y_train,X_l_train,X_nl_train,ind_l,ind_nl,U_train,X_l_pred,X_nl_pred,U_pred,B_pred,V,Tr,d,r,lambda,family,off, r_theta,eps)

    if (method == "Jackknife"){
      mfit=gplsvcm_est(Y_train,X_l_train,X_nl_train,U_train,V,Tr,d,r,lambda,family,off, r_theta,eps)
      if (length(ind_nl)==0){
        beta_hat=as.matrix(mfit$coefficients)
        alpha_pred=as.matrix(rep(0,dim(X_l_pred)[1]))
      }else{
        theta_hat=mfit$Qtheta
        alpha_pred=B_pred%*%theta_hat
      }
      if (length(ind_l)==0){
        beta_hat=c(0)
      }
      if (length(ind_nl)!=0 & length(ind_l)!=0){
        beta_hat=mfit$beta_hat
      }
      linkinv = family$linkinv
      ypred= linkinv(apply((X_nl_pred*alpha_pred),1,sum)+X_l_pred%*%beta_hat)     #mu_hat(X_test)
      q_res=quantile(lfit$residu,probs = cp,type=3)
      LL= ypred-q_res
      UU= ypred+q_res}

    if (method == "Jackknife+"){
      mat1=as.matrix(lfit$residu)
      mat2=mat1[, rep(1, dim(lfit$ypred)[2])]
      mat3=lfit$ypred-mat2
      mat4=lfit$ypred+mat2
      LL= apply(mat3,2,function(x) quantile(x,probs = 1-cp,type=3))
      UU= apply(mat4,2,function(x) quantile(x,probs = cp,type=3))}
  }
  if (family$family!="gaussian"){
    LL=ifelse(LL>=0,LL,0)
  }
  ###the prediction interval
  PIs=as.data.frame(cbind(LL,UU))
  colnames(PIs)=c("LL","UU")

  list(PIs=PIs)
}








