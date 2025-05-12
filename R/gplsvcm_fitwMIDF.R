#' Fitting the generalized partial linear spatially varying coefficient model
#' with Model Selection
#'
#' \code{gplsvcm_fitwMIDF} perform a model selection procedure to identify the
#' linear and nonlinear components first and then fit the corresponding
#' generalized partial linear spatially varying coefficient model.
#'
#' @importFrom MGLM kr
#' @importFrom BPST basis
#' @importFrom MASS glm.nb
#' @importFrom Matrix bdiag
#'
#' @param Y The response variable,a \code{n} by one matrix where \code{n} is the
#'   number of observations.
#' @param X The design matrix of \code{n} by \code{np} where \code{np} is the
#'   number of covariates (including intercept if intercept exists). Each row is a vector of the covariates for an
#'   observation.
#' @param U A \code{n} by two matrix where each row is the coordinates of an
#'   observation.
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
#' @param k_n The penalty parameter used in the model selection criteria. It
#'   need to be supplied only when the argument \code{method} is set to
#'   \code{NULL} --default is NULL.
#' @param off The offset -- default is 0.
#' @param r_theta The vector of the upper and lower bound of an interval to
#'   search for an additional parameter \code{theta} for negative binomial
#'   scenario -- default is c(2,8).
#' @param eps The error tolerance for the Pearson estimate of the scale
#'   parameter, which is as close as possible to 1, when estimating an
#'   additional parameter \code{theta} for negative binomial scenario -- default
#'   is 0.01.
#' @param method The type of model selection criteria, options are "AIC","BIC"
#'   and NULL which correspond to \code{k_n=2}, \code{k_n=log(n)} and
#'   \code{k_n=k_n} respectively -- default is "BIC".
#'
#' @return The function returns a list of fitted object information from S3
#'   class "gplsvcm", see the items of the list from \code{\link{gplsvcm_fit}}.
#'
#' @details The \code{gplsvcm_fitwMIDF} function is used to fit a generalized
#'   partial linear spatially varying coefficient model when the linear and
#'   nonlinear parts of the design matrix \code{X} are not known before
#'   analysis. The construction of the polynomial spline functions is via
#'   \code{\link[BPST]{basis}}. It first perform a model selection based on
#'   Generalized Information Criterion (GIC) and output the selected model by
#'   specifying the parameters \code{ind_l} and \code{ind_nl} of the function
#'   \code{gplsvcm_fit}. Then the selected model is fitted by the function
#'   \code{gplsvcm_fit}.
#'
#' @references Zhang et al.(2010) Regularization Parameter Selections via
#'   Generalized Information
#'   Criterion.\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2009.tm08013}
#'
#' @examples
#' # Population:
#' family=poisson()
#' ngrid = 0.02
#'
#' # Data generation:
#' pop = Datagenerator(family, ngrid)
#' N=nrow(pop)
#'
#' # Triangulations and setup:
#' Tr = Tr0; V = V0; n = 1000; d = 2; r = 1;
#'
#' # set up for smoothing parameters in the penalty term:
#' lambda_start=0.0001; lambda_end=10; nlambda=10;
#' lambda=exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))
#'
#' # Generate Sample:
#' ind_s=sample(N,n,replace=FALSE)
#' data=as.matrix(pop[ind_s,])
#' Y=data[,1]; X=data[,c(6:9)]; U=data[,c(10:11)];
#'
#' # True coefficents
#' alpha=data[,c(2:3)]; beta=data[,c(4:5)];
#'
#' # Fit the model with model selection based on AIC:
#' mfit1 = gplsvcm_fitwMIDF(Y, X, U, V, Tr, d , r , lambda,family,k_n=NULL,
#' method="AIC",off = 0,r_theta = c(2, 8), eps= 0.01)
#'
#' # Fit the model with model selection based on BIC:
#' mfit2 = gplsvcm_fitwMIDF(Y, X, U, V, Tr, d , r , lambda,family,k_n=NULL,
#' method="BIC",off = 0,r_theta = c(2, 8), eps= 0.01)
#'
#' # prediction intervals:
#' ind_l=mfit2$ind_l; ind_nl=mfit2$ind_nl;
#' set.seed(123)
#' PIs=compute_PIs(Y,X,ind_l,ind_nl,U,X,U,V,Tr,d,r,lambda,family,off = 0,
#' r_theta = c(2, 8), eps= 0.01,method="CV+", cp=0.95, nfold = 10)
#'
#' # prediction:
#' Y_hat = gplsvcm_predict(mfit2, X, U)
#'
#' # k-fold cross-validation:
#' set.seed(123)
#' MSPE = cv_gplsvcm(Y,X,ind_l,ind_nl,U,V,Tr,d,r,lambda,family,off = 0,r_theta =
#' c(2, 8), eps= 0.01,nfold=10)
#'
#' # plot the estimated coefficients
#' gplsvcm_plot(mfit2,gridnumber=100,display=c(1,1),xlab=c("u1","u1"),
#' ylab=c("u2","u2"),main=c(expression(paste("The Estimated Surface for","
#' ",hat(alpha)[1])),expression(paste("The Estimated Surface for","
#' ",hat(alpha)[2]))))
#'
#' @export


gplsvcm_fitwMIDF=function(Y,X,U,V,Tr,d=2,r=1,lambda= 10^seq(-6, 6, by = 0.5),
                          family,k_n=NULL,method="BIC",off = 0, r_theta = c(2, 8), eps= 0.01)
{
# structure selection
  np=ncol(X)
  ind_l=c()
  X_l=as.matrix(X[,ind_l])
  ind_nl=1:np
  X_nl=as.matrix(X[,ind_nl])
  mfit0 = suppressWarnings(gplsvcm_est(Y,X_l,X_nl,U,V,Tr,d,r,lambda,family,off,r_theta,eps))
  df=mfit0$df
  lam=mfit0$lambda_sel
  dev=dev_est(X_l,X_nl,mfit0,family)
  n=length(mfit0$Y)
  if (method == "AIC"){
    k_n=2
  }else if (method == "BIC"){
    k_n=log(n)
  }else{
    k_n=k_n
  }
  gic_old=dev+k_n*df
  gic_new=gic_old
  jnl=0
  while(gic_new<=gic_old & length(ind_nl)>0){
    gic_old=gic_new
    if(jnl!=0){
      ind_l=c(ind_l,ind_nl[jnl])
      ind_nl=ind_nl[-jnl]
    }
    nnl=length(ind_nl)
    if(nnl==0)
        break
    gic_all=c()
    for(inl in 1:nnl){
     tmp_l=c(ind_l,ind_nl[inl])
     tmp_nl=ind_nl[-inl]
     X_l=as.matrix(X[,tmp_l])
     X_nl=as.matrix(X[,tmp_nl])
     mfit=gplsvcm_est(Y,X_l,X_nl,U,V,Tr,d,r,lam,family,off,r_theta,eps)
     if(dim(X_nl)[2]==0){
      df=n-mfit$df.residual
     }else{
      df=mfit$df
     }
     dev=dev_est(X_l,X_nl,mfit,family)
     gic=dev+k_n*df
     gic_all=c(gic_all,gic)
   }
   jnl=which.min(gic_all)
   gic_new=gic_all[jnl]
  }
  # refit the model
  if(length(ind_l)!=0){
  ind_l=ind_l[order(ind_l)]}
  ind_nl=setdiff(1:np,ind_l)
  mfit=gplsvcm_fit(Y,X,ind_l,ind_nl,U,V,Tr,d=d,r=r,lambda=lambda,family,off=off,r_theta=r_theta,eps=eps)
  return(mfit)
}

