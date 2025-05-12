
#' @importFrom BPST basis
normalize_basis=function(uvpop=NULL,U,index=NULL,V,Tr,d ,r){
  # Preparation for spline approximation over triangulation
 if (!is.null(uvpop) & !is.null(index)){
         Basis_full = suppressWarnings(basis(V, Tr, d, r,uvpop))
         Q2 = Basis_full$Q2
         B = Basis_full$B
         BQ2 = B %*% Q2
         Bm=apply(BQ2,2,mean)
         Bm=Bm/Bm[1]
         Ba=t(do.call(cbind,(lapply(BQ2[,1],`*`,Bm))))
         Bc=BQ2-Ba
         Bc=Bc[index,-1]
  }else{
        Basis_full = basis(V, Tr, d, r,U)
        Q2 = Basis_full$Q2
        B = Basis_full$B
        BQ2 = B %*% Q2
        Bm=apply(BQ2,2,mean)
        Bm=Bm/Bm[1]
        Ba=t(do.call(cbind,(lapply(BQ2[,1],`*`,Bm))))
        Bc=BQ2-Ba
        Bc=Bc[,-1]
       }
  return(Bc)
}
