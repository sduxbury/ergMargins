#' Function to assess multicollinearity in ERGM using variance inflation factors
#' Measure based on Duxbury (2018) available in Sociological Methods and Research
#'
#' @param my.ergm is an ergm object
#'



vif.ergm<-function(my.ergm){


  if(class(my.ergm)%in%"btergm"){
    data_mat<-my.ergm@effects
    corr5<-stats::cor(data_mat[!rownames(data_mat)%in%"edges",
                        !colnames(data_mat)%in%"edges"]) ##omit edges term
    beta<-btergm::coef(my.ergm)
  }else{

      #get correlation matrix
    if(class(my.ergm)%in%"mlergm"){
      cor.mat<-stats::cov2cor(solve(my.ergm$information_matrix))
      beta<-my.ergm$theta

    }else{
      cor.mat<-stats::cov2cor(my.ergm$covar) #calculate correlation matrix
      beta<-stats::coef(my.ergm)

    }

    #omit edges, assign names to matrix
    rownames(cor.mat)<-colnames(cor.mat)<-names(beta)
    corr5<-cor.mat[!rownames(cor.mat)%in%"edges",
                   !colnames(cor.mat)%in%"edges"]
  }

  corr5<-corr5[!is.na(corr5[1:nrow(corr5)]),]
  corr5<-corr5[,which(!is.na(corr5[1,1:ncol(corr5)]))]

  VIFS<-matrix(0,nrow=1,ncol=ncol(corr5))

  for(i in 1:ncol(corr5)){

    gvec<-as.vector(corr5[-c(i),i]) ##create vector of correlations between covariate of interest and other covariates in the model
    tgvec<-t(gvec)
    xcor<-solve(corr5[-c(i),-c(i)]) ##create square matrix of correlations between covariates in the model other than the one of interest
    Rsq<-tgvec%*%xcor%*%gvec
    VIFS[1,i]<-1/(1-Rsq)
  }

  colnames(VIFS)<-names(beta[!names(beta)%in%"edges"])

  if(class(my.ergm)%in%"btergm"){

    warning("VIFS for bootstrap TERGM based on model matrix, not the covariance matrix of the estimator. Benchmarks used for MCMC ML estimation may not apply.")

   }else{
    message("Higher values indicate greater correlation.\nVIF > 20 is concerning, VIF > 100 indicates severe collinearity.")

  }
  VIFS
}



