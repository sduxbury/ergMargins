#' Function to conduct marginal mediation analysis when direct effects are interactions
#'
#'computes indirect effect using average marginal effects with bootstrapped standard errors
#' @param restricted.model is the model without the mediator variable(s)
#' @param full.model is the full model including all mediating variables
#' @param mediator is a mediating varaible or a character vector containing the names of the mediator variables if joint mediation is of interest
#' @param var1 is the main effect for the interaction.
#' @param var2 is the moderating effect.
#' @param inter is the interaction term
#' @param joint tells whether to assess mediation for the total effect of the interaction
#' @param inter tells whether to assess mediation for only the moderating effect, e.g., the differences between levels fo an interaction
#' At least joint or inter must be specified for the function to work


ergm.mod.mma<-function(restricted.model,full.model,var1, var2, inter,mediator,
                   at.2=NULL,joint=F,int.eff=F){


  if(joint==F & int.eff==F){
    stop("Please specify whether interested in assessing mediation for the total interaction effect, or only the interaction effect. ")
  }
  if(!is.null(at.2) & joint==F & length(at.2)==1){
    stop("Cannot calculate second differences for moderated effect when only 1 level of at.2 is specified. Set joint=T to calculate composite indirect effect.")
  }
  if(!is.null(at.2) & int.eff==T & length(at.2)==1){
    message("Cannot calculate second differences for moderated effect when only 1 level of at.2 is specified. Only composite indirect effects will be calculated.")
  }


  if(!is.na(pmatch("nodematch",inter))){
      indicator<-0
      if(network::is.directed(restricted.model$network)){
          if(!is.na(pmatch('nodeicov',var1))|!is.na(pmatch('nodeocov',var1))){
             indicator<-1}
        }
      if(!network::is.directed(restricted.model$network) & !is.na(pmatch("nodecov",var1))) {indicator<-1}
      if(indicator==1){
        ##matched nodal characteristics are not a product term, so compute marginal effects
        #for var 1 and var 2, then use results to compute marignal effect for interaction

        message("NOTE: when nodematch is specified with continuous main effects, it does not necessarily vary when the main effects change. It can therefore be treated as a main effect. Results returned are for ergm.mma.")
       return(ergm.mma(restricted.model,full.model,mediator=mediator,
                 direct.effect=inter))


      }
    }


  tot.AME<-ergm.AME(restricted.model,var1=var1,var2=var2,inter=inter,at.2=at.2,return.dydx=T,return.at.2 = T)
  at.2<-tot.AME[[2]]
  if(length(at.2)>1){
    tot.sec.diff<-tot.AME[[1]]$`Second differences`
  }
   tot.dydx<-tot.AME[[1]]$`Marginal effects`
  tot.AME<-tot.AME[[1]]$`Average Marginal effects`
  t.names<-rownames(tot.AME)

  p.AME<-ergm.AME(full.model,var1=var1,var2=var2,inter=inter,at.2=at.2,return.dydx=T)
  if(length(at.2)>1){
    p.sec.diff<-p.AME$`Second differences`
  }
  p.dydx<-p.AME$`Marginal effects`
  p.AME<-p.AME$`Average Marginal effects`

  rownames(tot.AME)<-paste("total effect:",t.names)
  rownames(p.AME)<-paste("partial effect:",rownames(p.AME))

  if(length(at.2)==1){

    ind.AME<-matrix(0,nrow=1,ncol=4)
    rownames(ind.AME)<-paste("indirect effect:",t.names)
    colnames(ind.AME)<-c("AME","Delta SE","Z","P")
    ind.AME[1,1]<-tot.AME[1,1]-p.AME[1,1]

      #compute covariance for all indirect effects
      cov.ame<-stats::cor(p.dydx[[1]][,1],tot.dydx[[1]][,1])
      cov.ame<-2*cov.ame*tot.AME[1,2]*p.AME[1,2]

    ind.AME[1,2]<-sqrt(tot.AME[,2]^2+p.AME[,2]^2-cov.ame)
    ind.AME[1,3]<-ind.AME[1,1]/ind.AME[1,2]
    ind.AME[1,4]<-2*stats::pnorm(-abs(ind.AME[1,3]))

    AME.out<-rbind(tot.AME,p.AME,ind.AME)
    prop.med<-1-(p.AME[1,1]/tot.AME[1,1])
    comment(AME.out)<-paste("Proportion mediated = ",prop.med)
    return(AME.out)
  }




    #compute total interaction effect
  if(joint==T){
  ###indirect effect

  mma.me<-tot.AME[,1]-p.AME[,1]
  cov.ame<-list()

  for(i in 1:length(p.dydx)){
    #compute covariance for all indirect effects
  cov.ame[[i]]<-stats::cor(p.dydx[[i]],tot.dydx[[i]])[1,1]
  cov.ame[[i]]<-2*cov.ame[[i]]*tot.AME[i,2]*p.AME[i,2]
  }
  cov.ame<-as.vector(cov.ame,mode="numeric")
  mma.se<-sqrt(tot.AME[,2]^2+p.AME[,2]^2-cov.ame)
  mma.z<-mma.me/mma.se
  p.mma<-2*stats::pnorm(-abs(mma.z))

  ind<-signif(cbind(mma.me,mma.se,mma.z,p.mma),digits=5)
  colnames(ind)<-colnames(tot.AME)

  if(length(mediator)>1) {
    mediator<-paste(mediator,collapse=", ")}

   rownames(ind)<-paste("indirect effect:",inter," == ",at.2,"->",mediator)

    proportion.mediated<-1-(p.AME[,1]/tot.AME[,1])
    ind<-cbind(ind,proportion.mediated)



    out1<-list(indirect.effects=ind,
               total.effects=tot.AME,
               partial.effects=p.AME)

    #compute summary statistics

    out2<-matrix(signif(c(mean(ind[,1]),mean(abs(ind[,3])))),ncol=2)
    out2<-cbind(out2,2*stats::pnorm(-abs(out2[,2])))
    colnames(out2)<-c("Mean second difference","Mean absolute Wald","P")

    out2<-list(summary.stats=out2,
               marginal.effects=out1)
    if(int.eff==F){
      out<-out2
    }

  }

    #compute moderator only
  if(int.eff==T){
    #difference in second difference
  third.diff<-tot.sec.diff[,1]-p.sec.diff[,1]

  diff.length<-length(p.dydx)-1
  tot.diff.list<-list()
  p.diff.list<-list()
  cov.list<-vector()
    #get correlation between second differences
    for(i in 1:diff.length){
      k<-i+1
      tot.diff.list[[i]]<-tot.dydx[[k]]-tot.dydx[[i]]
      p.diff.list[[i]]<-p.dydx[[k]]-p.dydx[[i]]
      cov.list[i]<-stats::cor(tot.diff.list[[i]],p.diff.list[[i]])[1,1]

    }
  cov.list<-cov.list*tot.sec.diff[,2]*p.sec.diff[,2]

  third.diff.se<-sqrt(tot.sec.diff[,2]^2+p.sec.diff[,2]^2-2*cov.list)
  third.diff.z<-third.diff/third.diff.se
  third.diff.p<-2*stats::pnorm(-abs(third.diff.z))

  if(length(third.diff)==1){#if there is only 1 third difference, no need for summary statistics
    third.diffs<-matrix(signif(c(third.diff,third.diff.se,third.diff.z,third.diff.p)),1,4,byrow=T)
    colnames(third.diffs)<-c("Third difference","Delta SE","Wald Z","P")

    out.list<-list(third.diffs=third.diffs,
                   total.second.diffs=tot.sec.diff,
                   partial.second.diffs=p.sec.diff,
                   total.marginal.effects=tot.AME,
                   partial.marginal.effects=p.AME)

    }else{#get summary statistics

      third.diffs<-cbind(third.diff,third.diff.se,third.diff.z,third.diff.p)
      colnames(third.diffs)<-c("Third difference","Delta SE","Wald Z","P")

      summary.output<-matrix(c(mean(third.diffs[,1]),mean(abs(third.diffs[,3])),0),ncol=3)
      summary.output[,3]<-2*stats::pnorm(-abs(summary.output[,2]))
      colnames(summary.output)<-c("Mean third diff","|Mean Wald Z|","P")

      out.list<-list(summary.output=summary.output,
                     third.diffs=third.diffs,
                     total.second.diffs=tot.sec.diff,
                     partial.second.diffs=p.sec.diff,
                     total.marginal.effects=tot.AME,
                     partial.marginal.effects=p.AME)


     }
   if(joint==F){
    out<-out.list
   }
  }
    #if interested in both interaction and joint effects, combine output
  if(joint==T & int.eff==T){
  out<-list(joint.effect=out2,
            moderator.effect=out.list)
  }


  return(out)
}
