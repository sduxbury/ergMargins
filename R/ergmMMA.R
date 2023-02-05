#' Function to conduct marginal mediation analysis
#'
#'computes indirect effect using average marginal effects with bootstrapped standard errors
#' @param restricted.model is the model without the mediator variable(s)
#' @param full.model is the full model including all mediating variables
#' @param mediator is a mediating varaible or a character vector containing the names of the mediator variables if joint mediation is of interest
#' @param direct.effect is the direct effect of interest


ergm.mma<-function(restricted.model,full.model,direct.effect,mediator,
                   at.controls=NULL,
                   control_vals=NULL){

  tot.AME<-ergm.AME(restricted.model,direct.effect,return.dydx=TRUE,
                    at.controls=at.controls,control_vals=control_vals)
  tot.dydx<-tot.AME$dydx
  tot.AME<-tot.AME$AME

  p.AME<-ergm.AME(full.model,direct.effect,return.dydx=TRUE,
                  at.controls=at.controls,control_vals=control_vals)
  p.dydx<-p.AME$dydx
  p.AME<-p.AME$AME

    #ensure equal dimensions. Typically only necessary in large networks.
  if(length(tot.dydx)!=length(p.dydx)){

    if(length(tot.dydx)<length(p.dydx)){
      p.dydx<-p.dydx[1:length(tot.dydx)]
    }else{
      tot.dydx<-tot.dydx[1:length(p.dydx)]
    }

  }

  rownames(tot.AME)<-paste("total effect:",rownames(tot.AME))
  rownames(p.AME)<-paste("partial effect:",rownames(p.AME))

    ###indirect effect

    mma.me<-tot.AME[1,1]-p.AME[1,1]
    cov.ame<-2*stats::cor(p.dydx,tot.dydx)*tot.AME[,2]*p.AME[,2]
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
