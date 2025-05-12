#' Produces coefficient function plots for a fitted generalized partial linear
#' spatially varying coefficient model.
#'
#' \code{gplsvcm_plot} produces the plots of the estimated coefficient functions
#' from a fitted \code{gplsvcm} object.
#'
#' @importFrom graphics points plot
#' @importFrom plot3D image2D
#'
#' @details This function used package \code{Triangulation} and \code{plot3D},
#'   see \code{\link[Triangulation]{TriPlot}} and \code{\link[plot3D]{image2D}}.
#'
#' @param mfit A fitted \code{gplsvcm} object returned from function
#'   \code{gplsvcm_fit} or \code{gplsvcm_fitwMIDF}.
#' @param gridnumber The number of grid points on one range for plots -- default
#'   is 100.
#' @param display If supplied then it is the vector for specifying how to
#'   display the estimated surfaces for the coefficient functions, used in
#'   \code{par(mfrow=)}.
#' @param xlab If supplied then is the vector of characters where each element
#'   is the x label for the estimated surface of one coefficient function.
#' @param ylab If supplied then is the vector of characters where each element
#'   is the y label for the estimated surface of one coefficient function.
#' @param main If supplied then is the vector of characters where each element
#'   is the title for the estimated surface of one coefficient function.
#' @param ... other graphics parameters to pass on to plotting commands. See
#'   details in \code{\link[plot3D]{image2D}}.
#' @return None
#'
#' @examples
#' # See an example of gplsvcm_fitwMIDF.
#'
#' @export

gplsvcm_plot = function(mfit, gridnumber = 100,display=NULL,xlab=NULL,ylab=NULL,main=NULL,...){
  while (!is.null(dev.list()))  dev.off()
  triplot = mTriPlot(mfit$V, mfit$Tr)

  np_nl = ncol(mfit$X_nl)
  if (np_nl == 0)
    stop("No coefficient functions in this model")

  S1b = cbind(min(mfit$V[,1]), max(mfit$V[,1]))
  S2b = cbind(min(mfit$V[,2]), max(mfit$V[,2]))

  # Generate dist.
  dist = max(S1b[2]-S1b[1], S2b[2]-S2b[1])/gridnumber
  uu = seq(S1b[1], S1b[2], dist);
  vv = seq(S2b[1], S2b[2], dist);
  n1 = length(uu); n2 = length(vv)
  u = rep(uu,n2); v = rep(vv,rep(n1,n2));
  uvpop = cbind(u,v)

  xran = (range(S1b)[2]-range(S1b)[1])/10
  yran = (range(S2b)[2]-range(S2b)[1])/10

  # Generate population basis for plot
  B_pred = suppressWarnings(basis(mfit$V,mfit$Tr,mfit$d,mfit$r,uvpop,FALSE,FALSE))
  Ind_all=B_pred$Ind
  if (!is.null(display)){
  par(mfrow = display)}
  for(ii in 1:np_nl){
    alpha=matrix(NA,n1*n2,1)
    alpha[Ind_all]=B_pred$Bi%*%mfit$Qtheta[,ii]
    image2D(uu,vv,z=matrix(alpha,n1,n2), xlim=c(min(S1b)-xran,max(S1b)+xran),
            ylim=c(min(S2b)-yran,max(S2b)+yran),xlab=xlab[ii],ylab=ylab[ii],main=main[ii],...)
  }

}

mTriPlot<-function(V, Tr, col=1, lwd=1){
  if(ncol(V)!=2){
    stop("Matrix of vectice should have 2 columns.")
  }
  Vmin=apply(V,2,min); Vmax=apply(V,2,max); Vplot=rbind(Vmin,Vmax);
  tri.plot=plot(x=Vplot,col="white",xlab="",ylab="",axes=F,main=expression("Triangulation Plot"))
  apply(Tr,1,tplot,V0=V,col=col,lwd=lwd)
  invisible()
}

tplot<-function(V0, Tr0, col=1, lwd=1){
  V1=V0[Tr0[1],]
  V2=V0[Tr0[2],]
  V3=V0[Tr0[3],]
  lines(x=rbind(V1,V2),col=col,lwd=lwd)
  lines(x=rbind(V2,V3),col=col,lwd=lwd)
  lines(x=rbind(V3,V1),col=col,lwd=lwd)
}
