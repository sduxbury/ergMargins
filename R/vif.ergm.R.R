#' Function to assess multicollinearity in ERGM using variance inflation factors
#' Measure based on Duxbury (2018) available in Sociological Methods and Research
#'
#' @param my.ergm is an ergm object
#'



vif.ergm<-function(my.ergm){

  cor.mat<-cov2cor(my.ergm$covar) #calculate correlation matrix
  corr5<-cor.mat[-c(1),-c(1)] ##omit edges term

  corr5<-corr5[!is.na(corr5[1:nrow(corr5)]),]
  corr5<-corr5[,which(!is.na(corr5[1,1:ncol(corr5)]))]

  VIFS<-matrix(0,nr=1,nc=ncol(corr5))

  for(i in 1:ncol(corr5)){

    gvec<-as.vector(corr5[-c(i),i]) ##create vector of correlations between covariate of interest and other covariates in the model
    tgvec<-t(gvec)
    xcor<-solve(corr5[-c(i),-c(i)]) ##create square matrix of correlations between covariates in the model other than the one of interest
    Rsq<-tgvec%*%xcor%*%gvec
    VIFS[1,i]<-1/(1-Rsq)
  }

    colnames(VIFS)<-names(my.ergm$coef[-c(1)])
 message("Higher values indicate greater correlation.\nVIF > 20 is concerning, VIF > 100 indicates severe collinearity.")
  VIFS
}



