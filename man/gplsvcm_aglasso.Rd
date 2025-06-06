% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gplsvcm_aglasso.R
\name{gplsvcm_aglasso}
\alias{gplsvcm_aglasso}
\title{Fitting the generalized partial linear spatially varying coefficient model
with Variable Selection and Model Structure Identification by (Adaptive) group lasso.}
\usage{
gplsvcm_aglasso(
  Y,
  X,
  uvpop = NULL,
  U,
  index = NULL,
  V,
  Tr,
  d = 2,
  r = 1,
  penalty = "agrLasso",
  lambda1 = exp(seq(log(1e-04), log(1), length.out = 50)),
  lambda2 = exp(seq(log(1e-04), log(1), length.out = 50)),
  family,
  criteria = "BIC",
  lambda = 10^seq(-6, 6, by = 0.5),
  off = 0
)
}
\arguments{
\item{Y}{The response variable,a \code{n} by one matrix where \code{n} is the
number of observations.}

\item{X}{The design matrix (without intercept) of \code{n} by \code{np} where \code{np} is the
number of covariates. Each row is a vector of the covariates for an
observation.}

\item{uvpop}{The coordinates of population grid points over the domain, default is \code{NULL}.}

\item{U}{A \code{n} by two matrix where each row is the coordinates of an
observation.}

\item{index}{The row indexes of the observed data points \code{U} in the population grid points
\code{uvpop}.}

\item{V}{A \code{N} by two matrix of vertices of a triangulation, where
\code{N} is the number of vertices and each row is the coordinates for a
vertex.}

\item{Tr}{A \code{n_Tr} by three triangulation matrix, where \code{n_Tr} is
the number of triangles in the triangulation and each row is the indices of
vertices in \code{V}.}

\item{d}{The degree of piecewise polynomials -- default is 2.}

\item{r}{The smoothness parameter and \code{r} \eqn{<} \code{d} -- default is
1.}

\item{penalty}{The shrinkage method for variable selection and model structure
identification, options are "agrLasso" for adaptive group lasso and "grLasso" for
group lasso -- default is "agrLasso".}

\item{lambda1}{The sequence of lambda values for group lasso,also need to be
specified for computing the weight if using adaptive group lasso --
default is a grid of 50 lambda values that ranges uniformly on the log scale
over 0.0001 to 1.}

\item{lambda2}{The sequence of lambda values used in the adaptive part for
adaptive group lasso -- default is a grid of 50 lambda values that ranges
uniformly on the log scale over 0.0001 to 1.}

\item{family}{The family object which specifies the distribution and link to
use (see \code{\link{glm}} and \code{\link{family}}).}

\item{criteria}{A character string specifying the criteria
of selecting lambda for (adaptive) group lasso. "BIC" is to use traditional
 Bayesian Information Criteria, "CV" is to use 10-fold cross validation --
 default is "BIC".}

\item{lambda}{The vector of the candidates of smoothing penalty parameter
-- default is grid points of 10 to the power of a sequence from -6 to 6 by 0.5.}

\item{off}{The offset -- default is 0.}
}
\value{
The function returns a list of fitted object information from S3
  class "gplsvcm", see the items of the list from \code{\link{gplsvcm_fit}}.
}
\description{
\code{gplsvcm_aglasso} perform variable selection and model structure identification
at the same time to select the covariates with the linear and nonlinear effects respectively
out of a large number of covariates using BIC or 10 fold cross validation and then fit the corresponding
generalized partial linear spatially varying coefficient model.
}
\details{
The \code{gplsvcm_aglasso} function is used to fit a generalized
  partial linear spatially varying coefficient model when there is a large
  number of covariates and the linear and nonlinear parts of the design
  matrix \code{X} are not known before analysis. The construction of the
  polynomial spline functions is via \code{\link[BPST]{basis}}. It first
  perform a variable selection and model structure identification through
  adaptive group lasso via \code{\link[grpreg]{grpreg}} or
  \code{\link[grpreg]{cv.grpreg}} and output the selected model by
  specifying the parameters \code{ind_l} and \code{ind_nl} of the function
  \code{gplsvcm_fit}. Then the selected model is fitted by the function
  \code{gplsvcm_fit}.
}
\examples{
# Population:
family=gaussian()
ngrid = 0.02

# Data generation:
pop = Datagenerator(family, ngrid)
N=nrow(pop)

# Triangulations and setup:
Tr = Tr0; V = V0; n = 1000; d = 2; r = 1;

# set up for smoothing parameters in the penalty term:
lambda_start=0.0001; lambda_end=10; nlambda=10;
lambda=exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))

# Generate Sample:
ind_s=sample(N,n,replace=FALSE)
data=as.matrix(pop[ind_s,])
Y=data[,1]; X=data[,c(6:9)]; U=data[,c(10:11)];

# True coefficents
alpha=data[,c(2:3)]; beta=data[,c(4:5)];

# Fit the model with model selection based on AIC:
mfit1 = gplsvcm_fitwMIDF(Y, X, U, V, Tr, d , r , lambda,family,k_n=NULL,
method="AIC",off = 0,r_theta = c(2, 8), eps= 0.01)

# Fit the model with model selection based on BIC:
mfit2 = gplsvcm_fitwMIDF(Y, X, U, V, Tr, d , r , lambda,family,k_n=NULL,
method="BIC",off = 0,r_theta = c(2, 8), eps= 0.01)

# prediction intervals:
ind_l=mfit2$ind_l; ind_nl=mfit2$ind_nl;
set.seed(123)
PIs=compute_PIs(Y,X,ind_l,ind_nl,U,X,U,V,Tr,d,r,lambda,family,off = 0,
r_theta = c(2, 8), eps= 0.01,method="CV+", cp=0.95, nfold = 10)

# prediction:
Y_hat = gplsvcm_predict(mfit2, X, U)

# k-fold cross-validation:
set.seed(123)
MSPE = cv_gplsvcm(Y,X,ind_l,ind_nl,U,V,Tr,d,r,lambda,family,off = 0,r_theta =
c(2, 8), eps= 0.01,nfold=10)

# plot the estimated coefficients
gplsvcm_plot(mfit2,gridnumber=100,display=c(1,1),xlab=c("u1","u1"),
ylab=c("u2","u2"),main=c(expression(paste("The Estimated Surface for","
",hat(alpha)[1])),expression(paste("The Estimated Surface for","
",hat(alpha)[2]))))

}
\references{
Wood, S., & Wood, M. S. (2015). Package ‘mgcv’. R package version, 1, 29.

Breheny P (2016).grpreg: Regularization Paths for Regression Models with
Grouped Covari-ates.Rpackage version 3.0-2,
URLhttps://CRAN.R-project.org/packages=grpreg.

Wang L, Wang G, Li X, Mu J, Yu S, Wang Y, Kim M, Wang J (2019).
BPST: Smoothing viaBivariate Spline over Triangulation.Rpackage version 1.0,
 URLhttps://GitHub.com/funstatpackages/BPST.
}
