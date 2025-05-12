#' Estimation for GPLSCVMs
#'
#' This is an internal function of \code{GPLSVCM} which is used in function
#' \code{gplsvcm_fit} and \code{gplsvcm_fitwMIDF}.
#'
#' @details The construction of the polynomial spline functions is via
#'   \code{\link[BPST]{basis}}.
#'
#' @param Y The response variable,a \code{n} by one matrix where \code{n} is the
#'   number of observations.
#' @param X_nl The matrix of nonlinear covariates for observations, of dimension
#'   \code{n} by \code{np_nl} where \code{n} is number of observations and
#'   \code{np_nl} is the number of nonlinear covarites.
#' @param X_l The matrix of linear covariates for observations, of dimension
#'   \code{n} by \code{np_l} where \code{n} is number of observations and
#'   \code{np_l} is the number of linear covarites.
#' @param U A \code{n} by two matrix where each row is the coordinates of an
#'   observation.
#' @param V A \code{N} by two matrix of vertices of a triangulation, where
#'   \code{N} is the number of vertices and each row is the coordinates for a
#'   vertex.
#' @param Tr A \code{n_Tr} by three triangulation matrix, where \code{n_Tr} is
#'   the number of triangles in the triangulation and each row is the indices of
#'   vertices in \code{V}.
#' @param d The degree of piecewise polynomials.
#' @param r The smoothness parameter and \code{r} \eqn{<} \code{d}.
#' @param lambda The vector of the candidates of penalty parameter.
#' @param family The family object which specifies the distribution and link to
#'   use (see \code{\link{glm}} and \code{\link{family}}).
#' @param off The offset.
#' @param r_theta The vector of the upper and lower bound of an interval to
#'   search for an additional parameter \code{theta} for negative binomial
#'   scenario.
#' @param eps The error tolerance for the Pearson estimate of the scale
#'   parameter, which is as close as possible to 1, when estimating an
#'   additional parameter \code{theta} for negative binomial scenario.


gplsvcm_est=function(Y,X_l,X_nl,U,V,Tr,d,r, lambda, family, off, r_theta, eps)
{
  if(dim(X_nl)[2] == 0 | all(X_nl == 0)){
    warning("There is no nonlinear coefficient,this will fit a GLM")
    if(family$family != "nb_bps"){
      return(glm(V1~.-1,data=as.data.frame(cbind(Y,X_l)),family = family))}
    if(family$family == "nb_bps"){
      return(glm.nb(V1~.-1,data=as.data.frame(cbind(Y,X_l))))}
  }else{
    if(dim(X_l)[2] == 0 | all(X_l == 0)){
      warning("There is no linear coefficient,this will fit a GSVCM")}
    if(all(X_l == 0)){
      X_l=X_nl[,c()]
    }
    if(family$family == "gaussian"){
      return(gplsvcm_est_gaus(Y,X_l,X_nl,U,V,Tr,d,r,lambda))
    }
    if(family$family == "poisson"){
      return(gplsvcm_est_pos(Y,X_l,X_nl,U,V,Tr,d,r,lambda, family, off, theta = 0))
    }
    if(family$family == "nb_bps"){
      return(gplsvcm_est_nb(Y,X_l,X_nl,U,V,Tr,d ,r,lambda, family, off, r_theta, eps))
    }
  }
}


gplsvcm_est_gaus=function(Y,X_l,X_nl,U,V,Tr,d ,r,lambda){

  # Preparation for spline approximation over triangulation
  Basis_full = suppressWarnings(basis(V, Tr, d, r, U))
  K = Basis_full$K
  Q2 = Basis_full$Q2
  B = Basis_full$B
  BQ2 = B %*% Q2
  P = t(Q2) %*% K %*% Q2
  ind_inside = Basis_full$Ind.inside
  tria_all = Basis_full$tria.all
  Y = Y[ind_inside]
  X_nl = as.matrix(X_nl[ind_inside, ])
  X_l=as.matrix(X_l[ind_inside, ])
  U = as.matrix(U[ind_inside, ])

  n_obs = length(Y)
  np_l=ncol(X_l)
  np_nl = ncol(X_nl)
  J = dim(BQ2)[2]
  # The design matrix
  Z = as.matrix(cbind(X_l,kr(X_nl, BQ2, byrow=TRUE)))
  # The modified block diagonal penalty matrix
  D = as.matrix(bdiag(diag(0,np_l),as.matrix(kronecker(diag(rep(1, np_nl)), P))))

  n_lam = length(lambda)
  est_all = matrix(rep(0, (np_nl*J+np_l) * n_lam), ncol = n_lam)
  gcv_all = rep(0, n_lam)
  df_all=c()

  # estimate the basis coefficients for each lambda.
  for(inl in 1:n_lam){
    temp = t(Z) %*% Z + lambda[inl] * D
    est_old = solve(temp, t(Z) %*% Y)      #Estimate the coefficients
    est_all[,inl] = as.matrix(est_old)

    # Compute the GCV criterion
    M_inv = solve(temp)
    S_lam = Z %*% tcrossprod(M_inv,Z)
    df = sum(diag(S_lam))
    df_all=c(df_all,df)
    Y_hat = Z %*% est_old
    gcv = n_obs * sum((Y-Y_hat)^2) / (n_obs - df)^2
    gcv_all[inl] = gcv
  }
  jnl = which.min(gcv_all)
  gcv = gcv_all[jnl]
  df = df_all[jnl]
  lambda_sel = lambda[jnl]
  est = est_all[,jnl]
  if (np_l==0){
    theta_hat = matrix(est, J, np_nl)
    beta_hat=NA
  }else{
    theta_hat = matrix(est[-(1:np_l)], J, np_nl)
    beta_hat=matrix(est[1:np_l],ncol=1)}
  alpha_hat = as.matrix(BQ2 %*% theta_hat)
  Qtheta= as.matrix(Q2 %*% theta_hat)
  list(Y=Y,X_l=X_l,X_nl=X_nl,U=U,alpha_hat = alpha_hat, beta_hat=beta_hat,Qtheta = Qtheta,lambda_sel = lambda_sel,gcv = gcv,df = df,B = B,
       Q2 = Q2, K = K,ind_inside = ind_inside, tria_all = tria_all)
}


gplsvcm_est_pos=function(Y,X_l,X_nl,U,V,Tr,d ,r ,lambda,family,off,theta){
  # Environmental Variables
  variance = family$variance
  linkinv = family$linkinv
  linkfun = family$linkfun
  mu_eta = family$mu.eta
  initialize = family$initialize

  # Preparation for spline approximation over triangulation
  Basis_full = suppressWarnings(basis(V, Tr, d, r, U))
  K = Basis_full$K
  Q2 = Basis_full$Q2
  B = Basis_full$B
  BQ2 = B %*% Q2
  P = t(Q2) %*% K %*% Q2
  ind_inside = Basis_full$Ind.inside
  tria_all = Basis_full$tria.all
  Y = Y[ind_inside]
  X_nl = as.matrix(X_nl[ind_inside, ])
  X_l=as.matrix(X_l[ind_inside, ])
  U = as.matrix(U[ind_inside, ])

  n_obs = length(Y)
  np_l=ncol(X_l)
  np_nl = ncol(X_nl)
  J = dim(BQ2)[2]
  # The design matrix
  Z = as.matrix(cbind(X_l,kr(X_nl, BQ2, byrow=TRUE)))
  # The modified block diagonal penalty matrix
  D = as.matrix(bdiag(diag(0,np_l),as.matrix(kronecker(diag(rep(1, np_nl)), P))))

  n_lam = length(lambda)
  est_all = matrix(rep(0, (np_nl*J+np_l) * n_lam), ncol = n_lam)
  gcv_all = rep(0, n_lam)
  df_all=c()

  # Perform the PIRLS algorithm for each lambda.
  for(inl in 1:n_lam){
    mu_0 = Y + 0.1
    eta_0 = linkfun(mu_0)
    if(theta != 0) V = variance(mu_0, theta) else V = variance(mu_0)
    G_inv = mu_eta(eta_0)

    Y_j = (eta_0 - off) + (Y - mu_0) / G_inv    #Y tilde
    W_j = (G_inv^2) / V     #W
    W_j = as.vector(W_j)
    temp1 = as.matrix(W_j * Z)       #WZ
    temp2 = as.vector(W_j * Y_j)       #WY
    temp3 = t(Z) %*% temp1 + lambda[inl] * D     #M(lambda)
    est_old = solve(temp3, t(Z) %*% temp2)      #Estimate the coefficients


    step = 0
    delta1 = 1
    delta2 = 1
    # Update the coefficient estimates until convergence
    while(delta1 > 1e-5 & sum(is.infinite(delta2)) == 0 & step <= 10){
      step = step + 1
      eta = Z %*% est_old + off
      mu = linkinv(eta)
      # The variance matrix V.
      if(theta != 0) V = variance(mu, theta) else V = variance(mu)

      G_inv = mu_eta(eta)
      Y_j = (eta - off) + (Y - mu) / G_inv
      W_j = as.vector((G_inv^2) / V)
      temp1 = as.matrix(W_j * Z)
      temp2 = as.vector(W_j * Y_j)
      temp3 = t(Z) %*% temp1 + lambda[inl] * D
      est_new = solve(temp3, t(Z) %*% temp2)

      # The distance between updated coefficient estimates and the old one.
      delta1 = sqrt(mean((est_new - est_old)^2))

      eta_new = Z %*% est_new
      delta2 = exp(eta_new)
      # Avoid infinite coefficient estimates
      if(sum(is.infinite(delta2)) == 0){
        est_old = est_new
      }
    }
    est_all[,inl] = as.matrix(est_old)

    # Compute the GCV criterion
    M_inv = solve(temp3)
    S_lam = Z %*% tcrossprod(M_inv,temp1)
    df = sum(diag(S_lam))
    df_all=c(df_all,df)
    Y_hat = Z %*% est_old
    gcv = n_obs * sum(W_j * (Y_j-Y_hat)^2) / (n_obs - df)^2
    gcv_all[inl] = gcv
  }
  jnl = which.min(gcv_all)
  gcv = gcv_all[jnl]
  df = df_all[jnl]
  lambda_sel = lambda[jnl]
  est = est_all[,jnl]
  if (np_l==0){
    theta_hat = matrix(est, J, np_nl)
    beta_hat=NA
  }else{
    theta_hat = matrix(est[-(1:np_l)], J, np_nl)
    beta_hat=matrix(est[1:np_l],ncol=1)}
  alpha_hat = as.matrix(BQ2 %*% theta_hat)
  Qtheta= as.matrix(Q2 %*% theta_hat)
  list(Y=Y,X_l=X_l,X_nl=X_nl,U=U,alpha_hat = alpha_hat, beta_hat=beta_hat,Qtheta = Qtheta,lambda_sel = lambda_sel,gcv = gcv,df = df,B = B,
       Q2 = Q2, K = K,ind_inside = ind_inside, tria_all = tria_all)
}


gplsvcm_est_nb =function(Y,X_l,X_nl,U,V,Tr,d,r ,lambda, family, off, r_theta, eps)
{
  # Environmental Variables
  variance = family$variance
  linkinv = family$linkinv
  linkfun = family$linkfun
  mu_eta = family$mu.eta
  initialize = family$initialize

  # Compute BQ2
  Basis_full = basis(V, Tr, d, r, U,FALSE,FALSE)
  B = Basis_full$B

  # The upper and lower bound for the theta parameter of neg-binomial family.
  lower = r_theta[1]
  upper = r_theta[2]

  # Compute the Pearson estimate of the scale parameter given the upper bound of the theta parameter
  result1 = gplsvcm_est_pos(Y,X_l,X_nl,U,V,Tr,d,r ,lambda, family, off, upper)
  if (dim(X_l)[2]==0){
    Y_hat = linkinv(kr(result1$X_nl, B, byrow = TRUE) %*% as.vector(result1$Qtheta))
  }else{
    Y_hat = linkinv(kr(result1$X_nl, B, byrow = TRUE) %*% as.vector(result1$Qtheta)+result1$X_l%*%as.vector(result1$beta_hat))}
  tt= sum((result1$Y - Y_hat)^2 / variance(Y_hat, upper))
  df_res = length(result1$Y) - result1$df
  scale_est_U=tt/df_res

  # Compute the Pearson estimate of the scale parameter given the lower bound of the theta parameter
  result2 = gplsvcm_est_pos(Y,X_l,X_nl,U,V,Tr,d,r,lambda, family, off, lower)
  if (dim(X_l)[2]==0){
    Y_hat = linkinv(kr(result2$X_nl, B, byrow = TRUE) %*% as.vector(result2$Qtheta))
  }else{
    Y_hat = linkinv(kr(result2$X_nl, B, byrow = TRUE) %*% as.vector(result2$Qtheta)+result2$X_l%*%as.vector(result2$beta_hat))}
  tt = sum((result2$Y - Y_hat)^2 / variance(Y_hat, lower))
  df_res = length(result2$Y) - result2$df
  scale_est_L=tt/df_res

  # Select the theta parameter such that the Pearson estimate of the scale parameter is close to 1
  scale_est_M = 2
  step = 1
  if ( (scale_est_L - 1) > 0 & (scale_est_U - 1) > 0 ) {
    scale_est_M = 1
    result = result2
    middle = lower
  }
  if ( (scale_est_L - 1) < 0 & (scale_est_U - 1) < 0 ) {
    scale_est_M = 1
    result = result1
    middle = upper
  }
  while(abs(scale_est_M - 1) > eps& step <= 10){
    middle = (upper + lower) / 2
    result = gplsvcm_est_pos(Y,X_l,X_nl,U,V,Tr,d,r, lambda, family, off, middle)
    if (dim(X_l)[2]==0){
      Y_hat = linkinv(kr(result$X_nl, B, byrow = TRUE) %*% as.vector(result$Qtheta))
    }else{
      Y_hat = linkinv(kr(result$X_nl, B, byrow = TRUE) %*% as.vector(result$Qtheta)+result$X_l%*%as.vector(result$beta_hat))}
    tt = sum((result$Y - Y_hat)^2 / variance(Y_hat, middle))
    df_res = length(result$Y) - result$df
    scale_est_M=tt/df_res
    if((scale_est_M - 1) * (scale_est_U - 1) < 0){
      lower = middle
      scale_est_L = scale_est_M
    } else {
      upper = middle
      scale_est_U = scale_est_M
    }
    step = step + 1
  }

  list(Y=result$Y,X_l=result$X_l,X_nl=result$X_nl,U=result$U,alpha_hat = result$alpha_hat, beta_hat=result$beta_hat,Qtheta = result$Qtheta,
       lambda_sel = result$lambda_sel, gcv = result$gcv, df = result$df, theta = middle,B = result$B, Q2 = result$Q2, K = result$K,ind_inside = result$ind_inside, tria_all = result$tria_all)
}

