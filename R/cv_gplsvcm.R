#' Calculate the K-fold Cross Validation Mean Square Prediction Error from a
#' fitted generalized partial linear spatially varying coefficient model.
#'
#' \code{cv_gplsvcm} implements K-fold cross-validation from a fitted
#' \code{gplsvcm} object, and returns the mean squared prediction error (MSPE).
#'
#' @details The construction of the polynomial spline functions is via
#'   \code{\link[BPST]{basis}}.
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
#' @param ind_l The vector of the indexes which indicate the columns of linear
#'   covariates in \code{X}.
#' @param ind_nl The vector of the indexes which indicate the columns of
#'   nonlinear covariates in \code{X}.
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
#' @param off The offset -- default is 0.
#' @param r_theta The vector of the upper and lower bound of an interval to
#'   search for an additional parameter \code{theta} for negative binomial
#'   scenario -- default is c(2,8).
#' @param eps The error tolerance for the Pearson estimate of the scale
#'   parameter, which is as close  to 1, when estimating an additional parameter
#'   \code{theta} for negative binomial scenario -- default is 0.01.
#' @param nfold The number of folds for cross validation -- default is 10.
#'
#' @return The mean square prediction error (MSPE).
#'
#' @examples
#' # See an example of gplsvcm_fitwMIDF.
#'
#' @export



cv_gplsvcm =function(Y, X,ind_l,ind_nl,U, V, Tr, d = 2, r = 1,
                     lambda = 10^seq(-6, 6, by = 0.5),family, off = 0,
                     r_theta = c(2, 8), eps= 0.01,nfold = 10)
{
  linkinv = family$linkinv
  Ball = suppressWarnings(basis(V, Tr, d, r, U))
  Q2 = Ball$Q2
  ind_inside = Ball$Ind.inside
  Y = as.matrix(Y[ind_inside])
  X_nl = as.matrix(X[ind_inside, ind_nl])
  X_l = as.matrix(X[ind_inside, ind_l])
  U = as.matrix(U[ind_inside, ])

  n = length(Y)
  sfold = round(n/nfold)
  Test = sample(1:n)
  cv_error = c()

length_record=0
  for(ii in 1:nfold){
    if(ii < nfold){
      Test_set = sort(Test[((ii - 1) * sfold + 1):(ii * sfold)])
    }
    if(ii == nfold){
      Test_set = sort(Test[((ii - 1) * sfold + 1):n])
    }
    Train_set = setdiff(1:n, Test_set)

    Y_test = as.matrix(Y[Test_set])
    Y_train =as.matrix(Y[Train_set])
    X_l_test = as.matrix(X_l[Test_set, ])
    X_l_train = as.matrix(X_l[Train_set, ])
    X_nl_test = as.matrix(X_nl[Test_set, ])
    X_nl_train = as.matrix(X_nl[Train_set, ])
    U_test = as.matrix(U[Test_set, ])
    U_train = as.matrix(U[Train_set, ])

    mfit_ii = gplsvcm_est(Y_train, X_l_train, X_nl_train, U_train, V, Tr, d, r, lambda, family, off, r_theta, eps)

    if (length(ind_nl)==0){
      n_test=dim(X_l_test)[1]
      beta_test=as.matrix(mfit_ii$coefficients)
      alpha_test=as.matrix(rep(0,n_test))
    }else{
      Ball_test = suppressWarnings(basis(V, Tr, d, r, U_test,FALSE,FALSE))
      alpha_test=Ball_test$B%*%mfit_ii$Qtheta
    }
    if (length(ind_l)==0){
      beta_test=c(0)
    }
    if (length(ind_nl)!=0 & length(ind_l)!=0){
     beta_test=mfit_ii$beta_hat
    }

    eta = apply((X_nl_test*alpha_test),1,sum)+X_l_test%*%beta_test
    Ypred_ii = linkinv(eta)


    pred_error = sum((Y_test - Ypred_ii)^2)
    cv_error = c(cv_error, pred_error)
    length_record=length_record+length(Y_test)
  }
  sspe=sum(cv_error)
  mspe=sspe/length_record
  return(mspe)
}
