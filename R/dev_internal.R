

##deviance function for gsvcm
dev1=function(mfit,family){
  alpha_hat=mfit$alpha_hat
  if(family$family == "gaussian"){
    return(sum((apply((mfit$X_nl*alpha_hat),1,sum)-mfit$Y)^2))
  }
  if(family$family == "poisson"){
    mu_hat=exp(apply((mfit$X_nl*alpha_hat),1,sum))
    return(2*sum(mfit$Y*log((mfit$Y+0.0001)/mu_hat)-mfit$Y+mu_hat))
  }
  if(family$family == "nb_bps"){
    mu_hat=exp(apply((mfit$X_nl*alpha_hat),1,sum))
    theta=mfit$theta
    return(2*sum(mfit$Y*log((mfit$Y+0.0001)/mu_hat)-(mfit$Y+theta)*log((1+mfit$Y/theta)/(1+mu_hat/theta))))
  }
}

##deviance function for glm
dev2=function(X_l,mfit,family)
{
  Y=as.vector(mfit$y)
  beta_hat=mfit$coefficients
  if(family$family == "gaussian"){
    return(sum((X_l%*%beta_hat-Y)^2))
  }
  if(family$family == "poisson"){
    mu_hat=exp(X_l%*%beta_hat)
    return(2*sum(Y*log((Y+0.0001)/mu_hat)-Y+mu_hat))
  }
  if(family$family == "nb_bps"){
    mu_hat=exp(X_l%*%beta_hat)
    theta=mfit$theta
    return(2*sum(Y*log((Y+0.0001)/mu_hat)-(Y+theta)*log((1+Y/theta)/(1+mu_hat/theta))))
  }
}

##deviance function for gplsvcm
dev3=function(X_l,X_nl,mfit,family)
{
  Y=as.vector(mfit$Y)
  alpha_hat=mfit$alpha_hat
  beta_hat=mfit$beta_hat
  if(family$family == "gaussian"){
    return(sum((apply((mfit$X_nl*alpha_hat),1,sum)+mfit$X_l%*%beta_hat-Y)^2))
  }
  if(family$family == "poisson"){
    mu_hat=exp(apply((mfit$X_nl*alpha_hat),1,sum)+mfit$X_l%*%beta_hat)
    return(2*sum(Y*log((Y+0.0001)/mu_hat)-Y+mu_hat))
  }
  if(family$family == "nb_bps"){
    mu_hat=exp(apply((mfit$X_nl*alpha_hat),1,sum)+mfit$X_l%*%beta_hat)
    theta=mfit$theta
    return(2*sum(Y*log((Y+0.0001)/mu_hat)-(Y+theta)*log((1+Y/theta)/(1+mu_hat/theta))))
  }
}

