#' Generating populations for simulation.
#'
#' \code{Datagenerator} is used to generate samples on horseshoe domain for
#' Scenario 1 (Gaussian), Scenario 2 (Poisson) and Scenario 3 (negative
#' binomial).
#'
#' @import mgcv
#' @importFrom MASS rnegbin
#' @importFrom boot inv.logit
#'
#' @param family The family object, specifying the distribution and link to use.
#'   Choose "gaussian" for Gaussian distribution, "poisson" for poisson
#'   distribution, and "nb_bps" for negative binomial distribution.
#' @param ngrid The distance between grid points -- default is set to 0.02.
#'
#' @return A data matrix with a response ('y'), true coefficient functions ('m1'
#'   and 'm2'), nonlinear covariates ('x1' and 'x2'), linear covariates ('x3'
#'   and 'x4') and locations ('u' and 'v').
#'
#' @details This function used the package \code{mgcv}, see
#'   \code{\link[mgcv]{fs.boundary}} and \code{\link[mgcv]{fs.test}}
#'
#' @examples
#' family=nb_bps()
#' ngrid = 0.02
#' pop = Datagenerator(family, ngrid)
#'
#' @export
#'
Datagenerator=function(family, ngrid){
  # Construct the HorseShoe boundary
  fsb=list(fs.boundary())
  # simulate some fitting data, inside boundary
  uu=seq(-1,3.5, ngrid)
  vv=seq(-1,1, ngrid)
  n1=length(uu)
  n2=length(vv)
  u=rep(uu,n2)
  v=rep(vv,rep(n1,n2))
  m1=fs.test(u,v,b=1)
  N=length(m1)
  m2=m1
  m2[!is.na(m1)]=mapply(beta1_geo, round(u[!is.na(m1)],2), round(v[!is.na(m1)],2))
  r1=-1
  r2=2

  x1=runif(N,-1,1)
  x2=runif(N,-1,1)
  x3=runif(N,-1,1)
  x4=runif(N,-1,1)


  # Scenario 1. Gaussian
  if(family$family=="gaussian"){
    sigma=0.5
    eps=suppressWarnings(rnorm(N,mean=0,sd=sigma))
    y=-2+m1*x1+m2*x2+r1*x3+r2*x4+eps

  }

  # Scenarios 2 and 3. Poisson and negative binomial
  if(family$family=="nb_bps"){
    m1=m1
    m2=m2
    log_mu=1+m1*x1+m2*x2+r1*x3+r2*x4
    mu=exp(log_mu)
    y=mu
    for(i in 1:N){
      y[i]=suppressWarnings(rnegbin(1,mu=mu[i],theta=6))
    }
  }
  if(family$family=="poisson"){
    m1=m1
    m2=m2
    log_mu=-2+m1*x1+m2*x2+r1*x3+r2*x4
    mu=exp(log_mu)
    y=mu
    for(i in 1:N){
      y[i]=suppressWarnings(rpois(1,mu[i]))
    }
  }

  if(family$family=="binomial"){
    m1=m1
    m2=m2
    logit_mu=1+m1*x1+m2*x2+r1*x3+r2*x4
    mu=inv.logit(logit_mu)
    y=mu
    for(i in 1:N){
      y[i]=rbinom(1,size=1,prob=mu[i])
    }
  }
  pop=cbind(y,m1,m2,r1,r2,x1,x2,x3,x4,u,v)
  pop=pop[!is.na(pop[,'m1']),]
  return(pop)
}


# modified beta 1 coefficient function:
beta1_geo=function(x,y){
  z=1.5*sin(pi/4*x*y)
  return(z)
}


