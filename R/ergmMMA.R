#' Function to conduct marginal mediation analysis
#'
#'computes indirect effect using average marginal effects with bootstrapped standard errors
#' @param restricted.model is the model without the mediator variable(s)
#' @param full.model is the full model including all mediating variables
#' @param mediator is a mediating varaible or a character vector containing the names of the mediator variables if joint mediation is of interest
#' @param direct.effect is the direct effect of interest


ergm.mma<-function(restricted.model,full.model,direct.effect,mediator,
                   at.controls=NULL,
                   control_vals=NULL,
                   ME="AME"){

  if(!ME%in%c("AME","MEM")){
    warning("ME must be specified as AME or MEM. Returning the AME.")
    ME<-"AME"
  }

    if(class(restricted.model)%in%"mlergm"){
      theta1<-restricted.model$theta
      theta2<-full.model$theta

    }else{
      theta1<-btergm::coef(restricted.model)
      theta2<-btergm::coef(full.model)

    }

  ##check at.controls appear in both models
  if(!is.null(at.controls)){
    if(!all(at.controls%in%names(theta1))){
      stop("Variables specified in any.controls must appear in both models to be used.")
    }
    if(!all(at.controls%in%names(theta2))){
      stop("Variables specified in any.controls must appear in both models to be used.")
    }
  }


  if(ME%in%"AME"){

    tot.AME<-ergm.AME(restricted.model,direct.effect,return.dydx=TRUE,
                    at.controls=at.controls,control_vals=control_vals)
    tot.Jac<-tot.AME$Jac
    tot.AME<-tot.AME$AME

    p.AME<-ergm.AME(full.model,direct.effect,return.dydx=TRUE,
                  at.controls=at.controls,control_vals=control_vals)
    p.Jac<-p.AME$Jac
    p.AME<-p.AME$AME
  }else{
    tot.AME<-ergm.MEM(restricted.model,direct.effect,return.dydx=TRUE,
                      at.controls=at.controls,control_vals=control_vals)
    tot.Jac<-tot.AME$Jac
    tot.AME<-tot.AME$MEM

    p.AME<-ergm.MEM(full.model,direct.effect,return.dydx=TRUE,
                    at.controls=at.controls,control_vals=control_vals)
    p.Jac<-p.AME$Jac
    p.AME<-p.AME$MEM

  }

    #ensure equal dimensions. Typically only necessary in large networks.
 # if(length(tot.dydx)!=length(p.dydx)){

  #  if(length(tot.dydx)<length(p.dydx)){
 #     p.dydx<-p.dydx[1:length(tot.dydx)]
  #  }else{
#      tot.dydx<-tot.dydx[1:length(p.dydx)]
  #  }

 # }


  ###estimate cross-model covariance
  names(p.Jac)<-names(theta2)
  names(tot.Jac)<-names(theta1)
  tot.Jac_expanded<-rep(0,length(p.Jac))
  names(tot.Jac_expanded)<-names(p.Jac)
  tot.Jac_expanded[names(tot.Jac_expanded)%in%names(tot.Jac)]<-tot.Jac

  #get vcovs
  if(class(restricted.model)%in%"mtergm"|class(restricted.model)%in%"btergm"){
    tot.vc <- stats::vcov(restricted.model@ergm)
    tot.vc <-tot.vc[!rownames(tot.vc)%in%"edgecov.offsmat",!colnames(tot.vc)%in%"edgecov.offsmat"]
    p.vc <- stats::vcov(full.model@ergm)
    p.vc <-p.vc[!rownames(p.vc)%in%"edgecov.offsmat",!colnames(p.vc)%in%"edgecov.offsmat"]

  }else{

    tot.vc <- stats::vcov(restricted.model)
    p.vc <- stats::vcov(full.model)

  }

  if(class(restricted.model)%in%"mlergm"){
    tot.vc<-solve(tot.vc)
    p.vc<-solve(p.vc)
  }


  tot_vc_expanded<-matrix(0,nrow=nrow(p.vc),ncol=ncol(p.vc))
  rownames(tot_vc_expanded)<-colnames(tot_vc_expanded)<-names(theta2)
  tot_vc_expanded[rownames(tot_vc_expanded)%in%names(theta1),
                colnames(tot_vc_expanded)%in%names(theta1)]<-tot.vc


  #calculate gradient and cross_mod covariance
  cross_mod_grad<-t(p.Jac)%*%t(as.matrix(tot.Jac_expanded))

  cross_mod_cov<-p.vc%*%cross_mod_grad%*%tot_vc_expanded


  rownames(tot.AME)<-paste("total effect:",rownames(tot.AME))
  rownames(p.AME)<-paste("partial effect:",rownames(p.AME))

    ###indirect effect

    mma.me<-tot.AME[1,1]-p.AME[1,1]
    cov.ame<-cross_mod_cov[direct.effect,direct.effect]
    mma.se<-sqrt(tot.AME[,2]^2+p.AME[,2]^2-cov.ame)
    mma.z<-mma.me/mma.se
    p.mma<-2*stats::pnorm(-abs(mma.z))

  ind<-matrix(signif(c(mma.me,mma.se,mma.z,p.mma),digits=5),nrow=1,ncol=4)
  colnames(ind)<-colnames(tot.AME)
   if(length(mediator)>1) {
    mediator<-paste(mediator,collapse=", ")}
  rownames(ind)<-paste("indirect effect:",direct.effect,"->",mediator)

  out<-rbind(tot.AME,p.AME,ind)
  proportion.mediated<-1-(p.AME[1]/tot.AME[1])
  attr(out,"description")<-paste("proportion of",direct.effect,"mediated by",mediator," = ",round(proportion.mediated,digits=3))


  return(out)
}
